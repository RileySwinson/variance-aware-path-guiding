#pragma once 

#if !defined(__MITSUBA_DS_UTIL_H_)
#define __MITSUBA_DS_UTIL_H_

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/random.h>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

MTS_NAMESPACE_BEGIN

struct Sample;
struct EnvironmentMap;
struct ErrorMetrics;

typedef boost::optional<EnvironmentMap> OptionalEnvMap;

struct MTS_EXPORT_CORE Sample {
    float value = 0.0f;
    
    float phi = 0.0f;
    float theta = 0.0f;
};

struct MTS_EXPORT_CORE EnvironmentMap {
	ref<Bitmap> bitmap;
	std::string filename;

	std::vector<std::pair<Point2i, Spectrum>> sampled_points;

	static OptionalEnvMap fetch(const boost::filesystem::path path)
	{
		if (path.extension() != ".exr" && path.extension() != ".hdr")
		{
			SLog(
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
					px[channel] -= noise; 
					if (px[channel] < 0.0f) px[channel] = 0.0f;
				}

				this->bitmap->setPixel(pt, px);
			}
		}
	}
};

struct MTS_EXPORT_CORE ErrorMetrics {
	static float MSE(Bitmap& bm1, Bitmap& bm2)
	{
		SAssert(bm1.getSize() == bm2.getSize());

		const int CHANNELS = 3;

		float err = 0.0f;
		for (int y = 0; y < bm1.getHeight(); ++y)
		{
			for (int x = 0; x < bm1.getWidth(); ++x)
			{
				Point2i pt(x, y);
				float bm1_rgb[CHANNELS];
				float bm2_rgb[CHANNELS];
				bm1.getPixel(pt).toLinearRGB(bm1_rgb[0], bm1_rgb[1], bm1_rgb[2]);
				bm2.getPixel(pt).toLinearRGB(bm2_rgb[0], bm2_rgb[1], bm2_rgb[2]);

				// Calculate error at given pixel by squaring the error
				float sum = 0.0f;
				for (int i = 0; i < CHANNELS; ++i)
				{
					sum += (bm1_rgb[i] - bm2_rgb[i]) * (bm1_rgb[i] - bm2_rgb[i]);
				}

				err += sum;
			}
		}

		// Divide by the pixel count
		err /= bm1.getPixelCount() * CHANNELS;
		return err;
	}

	static float RMSE(Bitmap& bm1, Bitmap& bm2)
	{
		// Take the square root of MSE for RMSE
		return std::sqrt(ErrorMetrics::MSE(bm1, bm2));
	}

	static float MAE(Bitmap& bm1, Bitmap& bm2)
	{
		SAssert(bm1.getSize() == bm2.getSize());

		const int CHANNELS = 3;

		float err = 0.0f;
		for (int y = 0; y < bm1.getHeight(); ++y)
		{
			for (int x = 0; x < bm1.getWidth(); ++x)
			{
				Point2i pt(x, y);
				float bm1_rgb[CHANNELS];
				float bm2_rgb[CHANNELS];
				bm1.getPixel(pt).toLinearRGB(bm1_rgb[0], bm1_rgb[1], bm1_rgb[2]);
				bm2.getPixel(pt).toLinearRGB(bm2_rgb[0], bm2_rgb[1], bm2_rgb[2]);

				// Calculate error at given pixel by taking the abs of the difference
				float sum = 0.0f;
				for (int i = 0; i < CHANNELS; ++i)
				{
					sum += std::abs(bm1_rgb[i] - bm2_rgb[i]);
				}

				err += sum;
			}
		}

		// Divide by the pixel count
		err /= bm1.getPixelCount() * CHANNELS;
		return err;
	}
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_DS_UTIL_H_ */