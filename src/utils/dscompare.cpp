#include <mitsuba/render/util.h>
#include <mitsuba/core/plugin.h>

#include <ds-compare/include.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/clamp.hpp>
#include <boost/functional/hash.hpp>

#include <unordered_map>

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
		cluster.attach(new VMFM());
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
			if (this->args.envmap_noise) envmap.noisify();
			envmap.precompute();

			/* Generate random samples and store them so they can be reused per data structure */
			std::vector<Sample> samples;
			uint32_t samples_learning = this->args.samples_learning;
			samples.reserve(samples_learning);
			for (uint32_t s_i = 0; s_i < samples_learning; ++s_i)
			{
				Point2 coords(random->nextFloat(), random->nextFloat());
				Sample sample = envmap.sample(this->args.mode, coords);
				if (this->args.samples_noise) sample.noisify();

				samples.push_back(sample);
			}

			/* Create folder for final output */
			const std::string base_name = this->args.result_path;
			const std::string folder_name = envmap.path.at(0);
			const std::string envmap_file_name = envmap.path.at(1);
			
			const std::string folder_path = base_name + "/" + folder_name + "/" + envmap_file_name;
			boost::filesystem::create_directories(folder_path);

			/* Generate bitmap for ground truth PDF */
			Float max = 0;
			EnvironmentMap gt = envmap
				.deep_copy(true)
				.map([&](Point2i coords, Spectrum& px) {
					Float lum = envmap.bitmap->getPixel(coords).getLuminance();
					for (int c = 0; c < envmap.bitmap->getChannelCount(); ++c) px[c] = lum;
					if (lum > max) max = lum;
				})
				.normalize(max);

			/* Write envmap to .exr file */
			gt.write(folder_path + "/base.exr");

			/* Initialize error metrics storage */
			std::vector<std::vector<float>> err_storage(cluster.largest() + 1);

			Log(EInfo, "Comparing data structures for envmap '%s'...", (folder_name + "/" + envmap_file_name).c_str());

			/* Iterate over data structures... */
			uint8_t curr_i = 1;
			cluster.for_each([&](DataStructure* ds) {
				if (ds->is_in(this->args.blacklist)) return;

				Log(EInfo, "[%s] | Active structure: '%s' (%i)",
					(std::to_string(curr_i) + "/" + std::to_string(cluster.size())).c_str(),
					ds->name().c_str(), ds->type()
				);

				/* Optional: Preprocess whatever has to be preprocessed per data structure */
				ds->preprocess();
				/* Store samples into the data structure */
				ds->store(samples);
				/* Optional: Postprocess whatever has to be postprocessed per data structure */
				ds->postprocess();

				/* Generate writable envmap with same properties as input envmap */
				Vector2i dims = envmap.bitmap->getSize();
				const std::string envmap_path = folder_path + "/" + std::to_string(ds->type()) + ".exr";

				Float max = 0;
				EnvironmentMap em = envmap
					.deep_copy(true)
					.map([&](Point2i coords, Spectrum& px) {
						Point2 norm((Float) coords.x / dims.x, (Float) coords.y / dims.y);
						Float density = ds->eval(norm);
						for (int c = 0; c < envmap.bitmap->getChannelCount(); ++c) px[c] = density;
						if (density > max) max = density;
					})
					.normalize(max);

				/* Write envmap to .exr file */
				em.write(envmap_path);

				/* Compute metrics and store them */
				err_storage.at(ds->type()) = std::vector<float>{
					ErrorMetrics::MSE(*gt.bitmap, *em.bitmap),
					ErrorMetrics::MAE(*gt.bitmap, *em.bitmap),
					ErrorMetrics::RMSE(*gt.bitmap, *em.bitmap)
				};

				/* Wipe data structure to clean state for further usage */
				ds->wipe();

				curr_i++;
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
		}

		cluster.clear();
		return EXIT_SUCCESS;
	}

	MTS_DECLARE_UTILITY()
private:
	DSArguments args;

	void handle_clargs(int argc, char** argv)
	{
		namespace p_opt = boost::program_options;

		try
		{
			BoostOptions desc("Options/Arguments");
			desc.add_options()
				// Utility
				("help,h", "Display help text.")
				// General
				("path,p", p_opt::value<std::string>(&this->args.path)->default_value("./data/tests/envmaps/"), "Path to envmap folder.")
				("result-path,rp", p_opt::value<std::string>(&this->args.result_path)->default_value("./data/results"), "Path to output folder.")
				("samples-learning,sl", p_opt::value<uint32_t>(&this->args.samples_learning)->default_value(1024), "Envmap sample count.")
				("samples-guiding,sg", p_opt::value<uint32_t>(&this->args.samples_guiding)->default_value(524288), "Reconstruction sample count.")
				("sample-mode,sm", p_opt::value<Sample::Mode>(&this->args.mode)->default_value(Sample::Mode::Cosine), "Envmap sampling mode.")
				("blacklist,b", p_opt::value<std::vector<int>>(&this->args.blacklist)->multitoken(), "List of data structure indices that won't be run.")
				// Noise
				("noisy-envmap,ne", p_opt::value<bool>(&this->args.envmap_noise)->default_value(false), "Noisify input envmap?")
				("noisy-samples,ns", p_opt::value<bool>(&this->args.samples_noise)->default_value(false), "Noisify learning samples?")
				// Spherical Harmonics
				("sh-bands,shb", p_opt::value<int>(&this->args.sh_bands)->default_value(5), "Number of Spherical Harmonic bands.")
				("sh-depth,shd", p_opt::value<int>(&this->args.sh_depth)->default_value(12), "Depth of Spherical Harmonics.")
				// DTree
				("dt-fracloss,dtl", p_opt::value<DTreeParams::EBsdfSamplingFractionLoss>(&this->args.dt_frac_loss)->default_value(DTreeParams::EBsdfSamplingFractionLoss::ENone), "Loss function during gradient descent.")
				("dt-dirfilter,dtf", p_opt::value<DTreeParams::EDirectionalFilter>(&this->args.dt_dir_filter)->default_value(DTreeParams::EDirectionalFilter::ENearest), "Directional filter for splatting radiance samples.")
				("dt-threshold,dtt", p_opt::value<Float>(&this->args.dt_threshold)->default_value(0.01), "Threshold for subdividing leaf nodes (percentage).")
				("dt-iter,dti", p_opt::value<int>(&this->args.dt_iterations)->default_value(-1), "Stop after nth iteration, starting at 0 (-1 to disable).")
				("dt-max-depth,dtd", p_opt::value<int>(&this->args.dt_max_depth)->default_value(20), "Maximum depth.")
				// TileCoding
				("tilings,t", p_opt::value<int>(&this->args.tilings)->default_value(4), "Number of tilings.")
				("tiles-x,tx", p_opt::value<int>(&this->args.tiles_x)->default_value(16), "Number of tiles in x direction.")
				("tiles-y,ty", p_opt::value<int>(&this->args.tiles_y)->default_value(16), "Number of tiles in y direction.")
				// VMM
				("vmm-components,vc", p_opt::value<uint32_t>(&this->args.vmf_components)->default_value(8), "Number of initial VMM components.");

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
