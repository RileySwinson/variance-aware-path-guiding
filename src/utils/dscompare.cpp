#include <mitsuba/render/util.h>
#include <mitsuba/core/plugin.h>

#include <mitsuba/ds/include.h>

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/clamp.hpp>

#include <vector>

MTS_NAMESPACE_BEGIN

struct Arguments {
	std::string path = "";
	bool noisify = false;
	int sh_bands = 3;
};

struct EnvironmentMap {
	ref<Bitmap> bitmap;
	std::string filename;

	/*Spectrum sample(ref<Random> random)
	{
		Point2f uv(random->nextFloat(), random->nextFloat());
		return this->bitmap->getPixel(uv);
	}*/

	/// Noisifies the underlying bitmap via a custom Gaussian noise implementation
	void noisify()
	{
		const float STD_DEV = 0.1f;
		const float MEAN 	= 0.0f;
		ref<Random> random = new Random();

		for (int y = 0; y < bitmap->getHeight(); ++y)
		{
			for (int x = 0; x < bitmap->getWidth(); ++x)
			{
				if (random->nextFloat() >= 0.1) continue; // only noisify ~10% of the pixels

				auto pt = Point2i(x, y);
				Spectrum px = bitmap->getPixel(pt);
				
				for (int channel = 0; channel < 3; ++channel)
				{
					float noise = random->nextFloat() * STD_DEV + MEAN;
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
		this->args.path = "./data/tests/envmaps";
		this->args.noisify = false;

		/* Register data structures */
		DSCluster& cluster = DSCluster::get();
		cluster.attach(new SphericalHarmonics());
		cluster.attach(new TileCoding());
		cluster.attach(new DirectionalTree());
		cluster.attach(new GMM());
		// ^^^ ... append your data structures here as you please

		/* Construct data structures as needed */
		cluster.for_each([](DataStructure* ds) {
			if (ds->type() == DSType::DS_Invalid) return;
			ds->construct();
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

			std::cout << envmap.filename << std::endl;

			cluster.for_each([](DataStructure* ds) {
				// TODO:
				// [ ] Importance sampling
				// [ ] Spherical harmonics
				// [ ] Generate envmap from that
				// [ ] RMSE for now
				// [ ] Store image

				
			});
		}

		return 0;
	}

	MTS_DECLARE_UTILITY()
private:
	Arguments args;

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
