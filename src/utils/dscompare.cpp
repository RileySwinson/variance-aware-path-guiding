#include <mitsuba/render/util.h>
#include <mitsuba/core/plugin.h>

#include <mitsuba/ds/include.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/clamp.hpp>
#include <boost/functional/hash.hpp>

#include <vector>
#include <unordered_map>
#include <exception>

MTS_NAMESPACE_BEGIN

typedef boost::program_options::options_description BoostOptions;
typedef boost::program_options::variables_map BoostOptionsMap;

class DSComparer : Utility {
public:
	int run(int argc, char** argv)
	{
		/* Deal with CL arguments */
		handle_clargs(argc, argv);

		/* Register data structures */
		DSCluster& cluster = DSCluster::get();
		cluster.attach(new SphericalHarmonics());
		cluster.attach(new TileCoding());
		cluster.attach(new DirectionalTree());
		cluster.attach(new GMM());
		// ^^^ ... append your data structures here as you please

		/* Construct data structures as needed */
		cluster.for_each([&](DataStructure* ds) {
			if (ds->type() == DSType::DS_Invalid) return;
			ds->construct(this->args);
		});

		/* Initialize random generator */
		ref<Random> random = new Random();

		/* Iterate over all environment maps */
		for (const auto& entry : boost::filesystem::recursive_directory_iterator(this->args.path))
		{
			if (boost::filesystem::is_directory(entry)) continue;

			OptionalEnvMap fetched_envmap = EnvironmentMap::fetch(entry.path());
			if (!fetched_envmap) continue;

			EnvironmentMap envmap = fetched_envmap.get();
			if (this->args.noisify)
			{
				envmap.noisify();
			}

			/* Optional: Preprocess whatever has to be preprocessed per data structure */
			cluster.for_each([&](DataStructure* ds) {
				if (ds->type() == DSType::DS_Invalid) return;
				ds->preprocess();
			});

			/* Generate N samples and store them into each data structure */
			for (uint32_t s_count = 0; s_count < this->args.samples; ++s_count)
			{
				Point2f rnd(random->nextFloat(), random->nextFloat());
				Sample sample = envmap.sample(rnd);

				cluster.for_each([&](DataStructure* ds) {
					if (ds->type() == DSType::DS_Invalid) return;
					ds->store(sample);
				});
			}

			/* Optional: Postprocess whatever has to be postprocessed per data structure */
			cluster.for_each([&](DataStructure* ds) {
				if (ds->type() == DSType::DS_Invalid) return;
				ds->postprocess();
			});

			ref<Bitmap> bm = new Bitmap(
				Bitmap::EPixelFormat::ERGB,
				Bitmap::EComponentFormat::EFloat32,
				envmap.bitmap->getSize(),
				3,
				nullptr
			);

			// Key: sample pos
			// Value: sample values
			std::unordered_map<
				std::pair<int, int>, 
				std::vector<float>, 
				boost::hash<std::pair<int, int>>
			> s_map;

			/* Sample approximated guiding distribution */
			for (uint32_t s_count = 0; s_count < this->args.samples; ++s_count)
			{
				Point2f rnd(random->nextFloat(), random->nextFloat());
				// Spherical harmonics for now
				auto sh = cluster.obtain(DSType::DS_SphericalHarmonics);
				Sample sample = sh->sample(rnd);

				// Normalize
				sample.phi *= INV_TWOPI;
				sample.theta *= INV_PI;

				Point2i pt(sample.phi * envmap.bitmap->getWidth(), sample.theta * envmap.bitmap->getHeight());
				s_map[std::pair<int, int>(pt.x, pt.y)].push_back(sample.value);
			}

			for (int y = 0; y < bm->getHeight(); ++y)
			{
				for (int x = 0; x < bm->getWidth(); ++x)
				{
					Point2i pt(x, y);
					auto entry = s_map.find(std::pair<int, int>(pt.x, pt.y));
					if (entry == s_map.end()) continue;

					float l = 0.0f;
					int sample_count = entry->second.size();
					for (int i = 0; i < sample_count; ++i)
					{
						float pdf = entry->second.at(i);
						l += pdf / sample_count;
					}

					Spectrum px = bm->getPixel(pt);
					Spectrum sampled_px = envmap.bitmap->getPixel(pt);
					
					px = sampled_px * l;
					bm->setPixel(pt, px);
				}
			}

			ref<Bitmap> bm2 = new Bitmap(Bitmap::EPixelFormat::ERGB, Bitmap::EComponentFormat::EFloat32, envmap.bitmap->getSize(), 3, nullptr);

			for (int y = 0; y < bm2->getHeight(); ++y)
			{
				for (int x = 0; x < bm2->getWidth(); ++x)
				{
					Point2i pt(x, y);
					Spectrum px = bm2->getPixel(pt);

					px[0] = 0; px[1] = 0; px[2] = 0;
					bm2->setPixel(pt, px);
				}
			}

			for (auto e_pair : envmap.sampled_points)
			{
				Point2i pt = e_pair.first;
				Spectrum sampled_px = envmap.bitmap->getPixel(pt);
				Spectrum px = bm2->getPixel(pt);

				px = sampled_px;

				bm2->setPixel(pt, px);
			}

			std::cout << "RMSE: " << ErrorMetrics::RMSE(*envmap.bitmap, *bm) << std::endl;
			std::cout << "RMSE: " << ErrorMetrics::RMSE(*envmap.bitmap, *bm2) << std::endl;

			envmap.bitmap->write(Bitmap::EFileFormat::EOpenEXR, "original.exr");
			bm->write(Bitmap::EFileFormat::EOpenEXR, "result.exr");
			bm2->write(Bitmap::EFileFormat::EOpenEXR, "result2.exr");
			break;
		}

		return 0;
	}

	MTS_DECLARE_UTILITY()
private:
	DSArguments args;

	void handle_clargs(int argc, char** argv)
	{
		try
		{
			BoostOptions desc("Options/Arguments");
			desc.add_options()
				("help,h", "Display help text.")
				("path,p", boost::program_options::value<std::string>(&this->args.path)->default_value("./data/tests/envmaps/"), "Path to envmap folder.")
				("samples,s", boost::program_options::value<uint32_t>(&this->args.samples)->default_value(8192), "Sample count.")
				("noisify,n", boost::program_options::value<bool>(&this->args.noisify)->default_value(false), "Noisify input envmap?")
				("sh-bands,b", boost::program_options::value<int>(&this->args.sh_bands)->default_value(3), "Number of Spherical Harmonic bands.");

			BoostOptionsMap op_map;
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), op_map);
			boost::program_options::notify(op_map);

			if (op_map.count("help"))
			{
				std::cout << desc << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		catch (std::exception& e)
		{
			std::cerr << e.what() << std::endl;
		}
	}
};

MTS_EXPORT_UTILITY(DSComparer, "Utility plugin for running comparison tests on a set of data structures")
MTS_NAMESPACE_END
