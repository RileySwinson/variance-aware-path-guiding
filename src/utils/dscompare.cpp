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

class DSComparer : public Utility {
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
		cluster.attach(new BinaryTileCoding());
		cluster.attach(new DirectionalTree());
		cluster.attach(new VMFM());
		// ^^^ ... append your data structures here as you please

		/* Construct data structures as needed */
		cluster.for_each([&](DataStructure* ds) {
			ds->construct(this->args);
		});

		/* Initialize random generator */
		ref<Random> random = new Random();

		/* Initialize error metrics storage */
		StatTrak& tracker = StatTrak::get();
		tracker.reserve(cluster.largest() + 1);

		/* Iterate over all environment maps */
		for (const auto& entry : boost::filesystem::recursive_directory_iterator(this->args.comparer.path))
		{
			if (boost::filesystem::is_directory(entry)) continue;

			/* Fetch environment map */
			OptionalEnvMap fetched_envmap = EnvironmentMap::fetch(entry.path());
			if (!fetched_envmap) continue;
			
			EnvironmentMap envmap = fetched_envmap.get();
			if (this->args.noise.envmap) envmap.noisify();
			envmap.precompute();

			/* Generate random samples and store them so they can be reused per data structure */
			std::vector<Sample> samples;
			uint32_t samples_learning = this->args.comparer.samples_learning;
			uint32_t samples_guiding = this->args.comparer.samples_guiding;
			samples.reserve(samples_learning);

			for (uint32_t s_i = 0; s_i < samples_learning; ++s_i)
			{
				Point2 coords(random->nextFloat(), random->nextFloat());
				Sample sample = envmap.sample(this->args.comparer.mode, coords);
				if (this->args.noise.samples) sample.noisify();

				samples.push_back(sample);
			}

			/* Create folder for final output */
			const std::string base_name = this->args.comparer.result_path;
			const std::string folder_name = envmap.path.at(0);
			const std::string envmap_file_name = envmap.path.at(1);
			
			const std::string folder_path = base_name + "/" + folder_name + "/" + envmap_file_name;
			boost::filesystem::create_directories(folder_path);

			/* Generate bitmap for ground truth PDF & write to file */
			Float max = 0;
			EnvironmentMap gt_map = envmap
				.deep_copy(true)
				.map([&](Point2i coords, Point3& px) {
					Float lum = envmap.get_pixel_luminance(coords);
					for (int c = 0; c < envmap.bitmap->getChannelCount(); ++c)
					{
						px[c] = lum;
					}

					if (lum > max) max = lum;
				});

			if (this->args.comparer.normalize)
			{
				gt_map.normalize(max);
			}

			Float gt_mean = gt_map.mean();
			gt_map.write(folder_path + "/gt.exr");

			Log(EInfo, "Comparing data structures for envmap '%s'...", (folder_name + "/" + envmap_file_name).c_str());

			/* Iterate over data structures... */
			uint8_t curr_i = 1;
			auto& blacklist = this->args.comparer.blacklist;
			cluster.for_each([&](DataStructure* ds) {
				if (ds->is_in(blacklist.begin(), blacklist.end())) return;

				std::string ds_counter = (std::to_string(curr_i) + "/" + std::to_string(cluster.size()));
				Log(EInfo, "[%s] | Active structure: '%s' (%i)",
					ds_counter.c_str(), ds->name().c_str(), ds->type()
				);

				tracker.follow(ds->type());
				tracker.timer_start("store");

				/* Optional: Preprocess whatever has to be preprocessed per data structure */
				ds->preprocess();
				/* Store samples into the data structure */
				ds->store(samples);
				/* Optional: Postprocess whatever has to be postprocessed per data structure */
				ds->postprocess();

				tracker.timer_end("store");

				/* Sample the base map using the approximation stored within the data structure and store the values for further MD calculation */
				std::vector<Float> observations;
				observations.reserve(samples_guiding);

				EnvironmentMap sample_map = envmap.deep_copy(true);
				for (int i = 0; i < samples_guiding; ++i)
				{
					Point2 rnd(random->nextFloat(), random->nextFloat());
					Sample sample = ds->sample(rnd);

					if (!sample.is_valid())
					{
						Log(EWarn, "%s Obtained invalid sample [φ: %f, θ: %f, p: %f] -- ignoring it!",
							(std::string(ds_counter.length() + 2, ' ') + " └").c_str(), sample.phi, sample.theta, sample.pdf
						);
						continue;
					}

					Point2 spherical(sample.phi, sample.theta);
					auto uv_coords = Converter::spherical_to_uv(spherical);
					auto im_coords = Converter::uv_to_image(uv_coords, envmap.bitmap->getSize());

					auto f_x = envmap.get_pixel_luminance(im_coords);
					auto p_x = sample.pdf;
					observations.push_back((f_x / p_x) * INV_FOURPI);

					Point3 col(std::max((Float) 0.0, 1 - 500 * sample.pdf), 0, std::min((Float) 1.0, 500 * sample.pdf));
					sample_map.set_pixel_rgb(im_coords, col);
				}
				sample_map.write(folder_path + "/" + std::to_string(ds->type()) + "_samples.exr");

				/* Evaluate function approximation per pixel and store the results in a new envmap */
				double d_sum = 0;
				EnvironmentMap eval_map = envmap
					.deep_copy(true)
					.map([&](Point2i coords, Point3& px) {
						Point2 uv = Converter::image_to_uv(coords, envmap.bitmap->getSize());
						double density = ds->eval(uv) * std::sin(Converter::uv_to_spherical(uv).y);
						
						for (int c = 0; c < envmap.bitmap->getChannelCount(); ++c)
						{
							px[c] = density;
						}

						d_sum += density;
					});

				d_sum = (d_sum / eval_map.bitmap->getPixelCount()) * (2 * M_PI * M_PI);
				if (d_sum < 0.99 || d_sum > 1.01)
				{
					Log(EWarn, "%s PDF does not properly integrate even within tolerable error margin... %f ∉ [0.99, 1.01]",
						(std::string(ds_counter.length() + 2, ' ') + " └").c_str(), d_sum
					);
				}

				/* Write envmap to .exr file */
				const std::string envmap_path = folder_path + "/" + std::to_string(ds->type()) + ".exr";
				eval_map.write(envmap_path);

				/* Compute metrics and store them */
				tracker.store(MD, ErrorMetrics::MD(observations, gt_mean));
				tracker.store(RMSE, ErrorMetrics::RMSE(gt_map, eval_map));
				tracker.store(MSE, ErrorMetrics::MSE(gt_map, eval_map));
				tracker.store(MAE, ErrorMetrics::MAE(gt_map, eval_map));

				tracker.store(Memory, ds->memory());

				/* Wipe data structure to clean state for further usage */
				ds->wipe();

				curr_i++;
			});

			/* Create metrics.csv */
			const std::string metrics_file_name = "metrics.csv";
			const std::string output_path = folder_path + "/" + metrics_file_name;

			tracker.write(output_path);
			tracker.reset();
		}

		tracker.clear();
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
				("path,p", p_opt::value<std::string>(&this->args.comparer.path), "Path to envmap folder.")
				("result-path,rp", p_opt::value<std::string>(&this->args.comparer.result_path), "Path to output folder.")
				("samples-learning,sl", p_opt::value<uint32_t>(&this->args.comparer.samples_learning), "Envmap sample count.")
				("samples-guiding,sg", p_opt::value<uint32_t>(&this->args.comparer.samples_guiding), "Reconstruction sample count.")
				("sample-mode,sm", p_opt::value<Sample::Mode>(&this->args.comparer.mode), "Envmap sampling mode.")
				("blacklist,b", p_opt::value<std::vector<int>>(&this->args.comparer.blacklist)->multitoken(), "List of data structure indices that won't be run.")
				("normalize,n", p_opt::value<bool>(&this->args.comparer.normalize), "Normalize?")
				// Noise
				("noisy-envmap,ne", p_opt::value<bool>(&this->args.noise.envmap), "Noisify input envmap?")
				("noisy-samples,ns", p_opt::value<bool>(&this->args.noise.samples), "Noisify learning samples?")
				// Spherical Harmonics
				("sh-bands,shb", p_opt::value<int>(&this->args.sh.bands), "Number of Spherical Harmonic bands.")
				("sh-depth,shd", p_opt::value<int>(&this->args.sh.depth), "Depth of Spherical Harmonics.")
				("sh-use-offset,sho", p_opt::value<bool>(&this->args.sh.use_offset), "Apply offset to SHs?")
				// DTree
				("dt-fracloss,dtl", p_opt::value<DTreeParams::EBsdfSamplingFractionLoss>(&this->args.dt.frac_loss), "Loss function during gradient descent.")
				("dt-dirfilter,dtf", p_opt::value<DTreeParams::EDirectionalFilter>(&this->args.dt.dir_filter), "Directional filter for splatting radiance samples.")
				("dt-threshold,dtt", p_opt::value<Float>(&this->args.dt.threshold), "Threshold for subdividing leaf nodes (percentage).")
				("dt-iter,dti", p_opt::value<int>(&this->args.dt.iterations), "Stop after nth iteration, starting at 0 (-1 to disable).")
				("dt-max-depth,dtd", p_opt::value<int>(&this->args.dt.max_depth), "Maximum tree depth.")
				// TileCoding
				("tilings,t", p_opt::value<int>(&this->args.tc.tilings), "Number of tilings.")
				("tiles-x,tx", p_opt::value<int>(&this->args.tc.tiles_x), "Number of tiles in x direction.")
				("tiles-y,ty", p_opt::value<int>(&this->args.tc.tiles_y), "Number of tiles in y direction.")
				// Binary Tile Coding
				("btc-tilings,bt", p_opt::value<int>(&this->args.btc.tilings), "Number of tilings.")
				("btc-tiles-x,btx", p_opt::value<int>(&this->args.btc.tiles_x), "Number of base tiles in x direction.")
				("btc-tiles-y,bty", p_opt::value<int>(&this->args.btc.tiles_y), "Number of base tiles in y direction.")
				("btc-max-depth,btd", p_opt::value<int>(&this->args.btc.max_depth), "Maximum tree depth.")
				("btc-threshold,btt", p_opt::value<Float>(&this->args.btc.subdiv_threshold), "Assumed population mean for T-Test.")
				// VMM
				("vmm-components,vc", p_opt::value<uint32_t>(&this->args.vmf.components), "Number of initial VMM components.")
				("vmm-mode,vr", p_opt::value<bool>(&this->args.vmf.use_ruppert), "Use ruppert implementation?");

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
