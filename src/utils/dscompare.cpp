#include <mitsuba/render/util.h>
#include <mitsuba/core/plugin.h>

#include <mitsuba/ds/include.h>

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/clamp.hpp>

#include <vector>

MTS_NAMESPACE_BEGIN

struct EnvironmentMap {
	ref<Bitmap> bitmap;
	std::string filename;

	std::vector<std::pair<Point2i, Spectrum>> sampled_points;

	Sample sample(Point2f& sample)
	{
		Vector dir = warp::squareToCosineHemisphere(sample);

		/* Transform to (hemi)spherical coordinates and normalize */
		float theta = std::acos(dir.z);
		float phi = std::atan2(dir.y, dir.x);
		if (phi < 0) phi += 2 * M_PI;

		Point2f uv_norm(
			phi * INV_TWOPI, 	// normalize phi into range [0.0, 1.0)
			theta * INV_PI 		// normalize theta into range [0.0, 0.5)
			// alt. for y: (1.0f - dir.z) * 0.5f // invert and halve (value is already in range [0.0, 1.0))
		);

		Point2i uv(
			uv_norm.x * this->bitmap->getWidth(),
			uv_norm.y * this->bitmap->getHeight()
		);

		/* UNCOMMENT FOR SAMPLE VISUALIZATION */
		sampled_points.push_back(std::pair<Point2i, Spectrum>(uv, this->bitmap->getPixel(uv)));

		Sample sample_data = {
			.value = this->bitmap->getPixel(uv).getLuminance(),
			.phi = phi,
			.theta = theta
		};

		/* Return found texel */
		return sample_data;
	}

	/// Noisifies the underlying bitmap via a custom Gaussian noise implementation
	void noisify(float noise_perc = 0.2f)
	{
		const float STD_DEV = 0.1f;
		const float MEAN 	= 0.0f;
		ref<Random> random = new Random();

		for (int y = 0; y < bitmap->getHeight(); ++y)
		{
			for (int x = 0; x < bitmap->getWidth(); ++x)
			{
				if (random->nextFloat() >= noise_perc) continue; // only alter a pixel with a certain chance

				auto pt = Point2i(x, y);
				Spectrum px = bitmap->getPixel(pt);
				
				float noise = random->nextFloat() * STD_DEV + MEAN;

				for (int channel = 0; channel < 3; ++channel)
				{
					px[channel] = ((px[channel] - noise) < 0.0f) ? 0.0f : (px[channel] - noise); // subtract, we want to have slightly darker pixels
				}

				this->bitmap->setPixel(pt, px);
			}
		}
	}
};

typedef boost::optional<EnvironmentMap> OptionalEnvMap;

class DSComparer : Utility {
public:
	int run(int argc, char** argv)
	{
		/* Deal with CL arguments */
		// TODO: use boost program options or the build in fetching stuff (see kdbench.cpp)
		//this->args.path = boost::optional<std::string>(argv[1]).get_value_or("./data/tests/envmaps");
		this->args.path = "./data/tests/envmaps/hdrihaven";
		this->args.noisify = true;

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

			OptionalEnvMap fetched_envmap = fetch_envmap(entry.path());
			if (!fetched_envmap) continue;

			EnvironmentMap envmap = fetched_envmap.get();
			if (this->args.noisify)
			{
				envmap.noisify();
			}

			/* Optional: Preprocess whatever has to be preprocessed per data structure */
			cluster.for_each([&](DataStructure* ds) {
				ds->preprocess();
			});

			/* Generate N samples and store them into each data structure */
			for (uint32_t s_count = 0; s_count < this->args.samples; ++s_count)
			{
				Point2f rnd(random->nextFloat(), random->nextFloat());
				Sample sample = envmap.sample(rnd);

				cluster.for_each([&](DataStructure* ds) {
					ds->store(sample);
				});
			}

			/* Optional: Postprocess whatever has to be postprocessed per data structure */
			cluster.for_each([&](DataStructure* ds) {
				ds->postprocess();
			});

			
			for (int y = 0; y < envmap.bitmap->getHeight(); ++y)
			{
				for (int x = 0; x < envmap.bitmap->getWidth(); ++x)
				{
					Point2i pt(x, y);
					Spectrum px = envmap.bitmap->getPixel(pt);

					px[0] = 0; px[1] = 0; px[2] = 0;
					envmap.bitmap->setPixel(pt, px);
				}
			}

			/* Sample approximated guiding distribution */
			for (uint32_t s_count = 0; s_count < this->args.samples; ++s_count)
			{
				Point2f rnd(random->nextFloat(), random->nextFloat());
				// Spherical harmonics for now
				auto sh = cluster.obtain(DSType::DS_SphericalHarmonics);
				Sample sample = sh->sample(rnd);

				//std::cout << sample.phi << " " << sample.theta << std::endl;

				// Normalize
				sample.phi *= INV_TWOPI;
				sample.theta *= INV_PI;

				Point2i pt(sample.phi * envmap.bitmap->getWidth(), sample.theta * envmap.bitmap->getHeight());
				Spectrum px = envmap.bitmap->getPixel(pt);

				px[0] = sample.value * 255;
				px[1] = sample.value * 255;
				px[2] = sample.value * 255;

				envmap.bitmap->setPixel(pt, px);
			}

			// TODO:
			// [X] Importance sampling
			// [ ] Spherical harmonics
			// [ ] Generate envmap from that
			// [ ] RMSE for now
			// [ ] Store image

			/*for (int y = 0; y < envmap.bitmap->getHeight(); ++y)
			{
				for (int x = 0; x < envmap.bitmap->getWidth(); ++x)
				{
					Point2i pt(x, y);
					Spectrum px = envmap.bitmap->getPixel(pt);

					px[0] = 0; px[1] = 0; px[2] = 0;
					envmap.bitmap->setPixel(pt, px);
				}
			}

			for (auto e_pair : envmap.sampled_points)
			{
				Point2i pt = e_pair.first;
				Spectrum px = envmap.bitmap->getPixel(pt);

				px[0] = e_pair.second.getLuminance();
				px[1] = e_pair.second.getLuminance();
				px[2] = e_pair.second.getLuminance();

				envmap.bitmap->setPixel(pt, px);
			}

			envmap.bitmap->write(Bitmap::EFileFormat::EOpenEXR, "test.exr");
			*/

			envmap.bitmap->write(Bitmap::EFileFormat::EOpenEXR, "test3.exr");
			break;
		}

		return 0;
	}

	MTS_DECLARE_UTILITY()
private:
	DSArguments args;

	/**
	 * 
	 * 
	 * @param path
	 * @return 
	 */
	OptionalEnvMap fetch_envmap(const boost::filesystem::path path)
	{
		if (path.extension() != ".exr" && path.extension() != ".hdr")
		{
			Log(
				ELogLevel::EWarn, "Environment map \"%s\" has invalid file format. Expected \".exr\" or \".hdr\", found: \"%s\" -- skipping.",
				path.string().c_str(), path.extension().c_str()
			);
			return boost::none;
		}

		EnvironmentMap envmap = {
			.bitmap = new Bitmap(path.string()),
			.filename = path.string()
		};

		return envmap;
	}
};

MTS_EXPORT_UTILITY(DSComparer, "Utility plugin for running comparison tests on a set of data structures")
MTS_NAMESPACE_END
