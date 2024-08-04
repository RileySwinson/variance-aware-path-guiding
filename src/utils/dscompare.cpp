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
		cluster.attach(new Unidirectional());
		cluster.attach(new SphericalHarmonics());
		cluster.attach(new TileCoding());
		cluster.attach(new DirectionalTree());
		cluster.attach(new GMM());
		// ^^^ ... append your data structures here as you please

		/* Construct data structures as needed */
		cluster.for_each([&](DataStructure* ds) {
			ds->construct(this->args);
		});

		/* Initialize random generator */
		ref<Random> random = new Random();

		/* Iterate over all environment maps */
		for (const auto& entry : boost::filesystem::recursive_directory_iterator(this->args.path))
		{
			if (boost::filesystem::is_directory(entry)) continue;

			/* Fetch environment map */
			OptionalEnvMap fetched_envmap = EnvironmentMap::fetch(entry.path());
			if (!fetched_envmap) continue;
			
			EnvironmentMap envmap = fetched_envmap.get();
			if (this->args.noisify) envmap.noisify();
			envmap.precompute();

			/* Generate random samples and store them so they can be reused per data structure */
			std::vector<Sample> samples;
			uint32_t samples_learning = this->args.samples_learning;
			for (uint32_t s_i = 0; s_i < samples_learning; ++s_i)
			{
				Point2 coords(random->nextFloat(), random->nextFloat());
				Sample sample = envmap.sample(this->args.mode, coords);

				samples.push_back(sample);
			}

			/* Create folder for final output */
			const std::string base_name = "./data/results";
			const std::string folder_name = envmap.path.at(0);
			const std::string envmap_file_name = envmap.path.at(1);
			
			const std::string folder_path = base_name + "/" + folder_name + "/" + envmap_file_name;
			boost::filesystem::create_directories(folder_path);

			/* Initialize error metrics storage */
			std::vector<std::vector<float>> err_storage(cluster.largest() + 1);

			Log(EInfo, "Comparing data structures for envmap '%s'...", (envmap.path.at(0) + "/" + envmap.path.at(1)).c_str());

			/* Iterate over data structures... */
			cluster.for_each([&](DataStructure* ds) {
				/* Optional: Preprocess whatever has to be preprocessed per data structure */
				ds->preprocess();
				/* Store samples into the data structure */
				ds->store(samples);
				/* Optional: Postprocess whatever has to be postprocessed per data structure */
				ds->postprocess();

				/* Samples will be stored in a map with key = pos, value = sample vector */
				std::unordered_map<
					std::pair<int, int>, 
					std::vector<Sample>, 
					boost::hash<std::pair<int, int>>
				> s_map;

				/* Sample approximated guiding distribution */
				uint32_t samples_guiding = this->args.samples_guiding;
				for (uint32_t s_i = 0; s_i < samples_guiding; ++s_i)
				{
					Point2 rnd(random->nextFloat(), random->nextFloat());
					Sample sample = ds->sample(rnd);

					// Normalize
					float u = sample.phi * INV_TWOPI;
					float v = sample.theta * INV_PI;

					Point2i pt(u * envmap.bitmap->getWidth(), v * envmap.bitmap->getHeight());
					//sample.value = envmap.bitmap->getPixel(pt).getLuminance();
					s_map[std::pair<int, int>(pt.x, pt.y)].push_back(sample);
				}

				/* Generate writable envmap with same properties as input envmap */
				ref<Bitmap> bm = envmap.gen_empty_bitmap();

				/* Fill envmap */
				for (int y = 0; y < bm->getHeight(); ++y)
				{
					for (int x = 0; x < bm->getWidth(); ++x)
					{
						Point2i pt(x, y);
						auto entry = s_map.find(std::pair<int, int>(pt.x, pt.y));
						if (entry == s_map.end()) continue;

						Spectrum px = bm->getPixel(pt);
						Spectrum sampled_px = envmap.bitmap->getPixel(pt);
						std::vector<Sample> samples = entry->second;

						Float base_value = 0.9;
						Float noise_factor = base_value + ((1 - base_value) * random->nextFloat());

						for (int channel = 0; channel < envmap.bitmap->getChannelCount(); ++channel)
						{
							Float value = sampled_px[channel];
							value *= noise_factor;
							if (value < 0) value = 0;
							
							Float l = 0;
							for (const auto& sample : samples)
							{
								l += value * (1 / sample.pdf); // f(x) * w(x) ... w(x) = 1 / pdf
							}

							l /= samples.size(); // (1 / N) * sum(...)
							px[channel] = l;
						}

						bm->setPixel(pt, px);
					}
				}

				/* Compute metrics and store them */
				err_storage.at(ds->type()) = std::vector<float>{
					ErrorMetrics::MSE(*envmap.bitmap, *bm),
					ErrorMetrics::MAE(*envmap.bitmap, *bm),
					ErrorMetrics::RMSE(*envmap.bitmap, *bm)
				};

				/* Write envmap bitmap to .exr file */
				const std::string envmap_path = folder_path + "/" + std::to_string(ds->type()) + ".exr";
				bm->write(Bitmap::EFileFormat::EOpenEXR, envmap_path);

				/* Wipe data structure to clean state for further usage */
				ds->wipe();
			});

			/* Create metrics.csv and fill it */
			std::ofstream output;
			const std::string metrics_file_name = "metrics.csv";
			const std::string output_path = folder_path + "/" + metrics_file_name;
			output.open(output_path, std::ios::out);
			
			output << ",MSE,RMSE,MAE,\n";

			std::string res;
			for (int i = 0; i < err_storage.size(); ++i)
			{
				res += "" + std::to_string(i) + ",";

				const auto metrics = err_storage.at(i);
				if (metrics.empty())
					res += ",,,";
				else
					for (const auto metric : metrics)
						res += std::to_string(metric) + ",";

				res += "\n";
			}

			output << res;

			/*for (int y = 0; y < envmap.bitmap->getHeight(); y += 2)
			{
				for (int x = 0; x < envmap.bitmap->getWidth(); x += 2)
				{
					Point2f rnd(
						((float) y / envmap.bitmap->getHeight()) * M_PI,
						((float) x / envmap.bitmap->getWidth()) * (2 * M_PI)
					);

					Point2f rnd2((float) x / envmap.bitmap->getWidth(), (float) y / envmap.bitmap->getHeight());

					// Spherical harmonics for now
					auto ds = cluster.obtain(DSType::DS_DTree);
					SphericalHarmonics* sh = dynamic_cast<SphericalHarmonics*>(ds);

					float result = sh->eval(rnd.x, rnd.y);

					Point2i pt(x, y);
					s_map[std::pair<int, int>(pt.x, pt.y)].push_back(result);
				}
			}*/
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
				// Utility
				("help,h", "Display help text.")
				// General
				("path,p", boost::program_options::value<std::string>(&this->args.path)->default_value("./data/tests/envmaps/"), "Path to envmap folder.")
				("samples-learning,sl", boost::program_options::value<uint32_t>(&this->args.samples_learning)->default_value(1024), "Envmap sample count.")
				("samples-guiding,sg", boost::program_options::value<uint32_t>(&this->args.samples_guiding)->default_value(524288), "Reconstruction sample count.")
				("sample-mode,sm", boost::program_options::value<Sample::Mode>(&this->args.mode)->default_value(Sample::Mode::Cosine), "Envmap sampling mode.")
				("noisify,n", boost::program_options::value<bool>(&this->args.noisify)->default_value(false), "Noisify input envmap?")
				// Spherical Harmonics
				("sh-bands,shb", boost::program_options::value<int>(&this->args.sh_bands)->default_value(5), "Number of Spherical Harmonic bands.")
				("sh-depth,shd", boost::program_options::value<int>(&this->args.sh_depth)->default_value(12), "Depth of Spherical Harmonics.");

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
