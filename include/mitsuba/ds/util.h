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
    Float value = 0.0f;
    Float phi = 0.0f;
    Float theta = 0.0f;

	enum Mode {
		Cosine,
		Native,
		Sphere
	};

	friend std::istream& operator>>(std::istream& in, Sample::Mode& mode)
	{
		std::string token;
		in >> token;

		if (token == "cosine")
		{
			mode = Sample::Mode::Cosine;
			return in;
		}
		if (token == "native")
		{
			mode = Sample::Mode::Native;
			return in;
		}
		if (token == "sphere")
		{
			mode = Sample::Mode::Sphere;
			return in;
		}

		in.setstate(std::ios_base::failbit);
		return in;
	};
};

struct MTS_EXPORT_CORE EnvironmentMap {
	ref<Bitmap> bitmap;
	std::string filename;
	bool precomputed = false;

	std::vector<std::pair<Point2i, Spectrum>> sampled_points;

	Float bitmap_integral;
	std::vector<Float> row_avgs;

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
			.filename = path.filename().string()
		};

		return envmap;
	}

	Sample sample(Sample::Mode mode, Point2f& sample)
	{
		switch(mode)
		{
			case Sample::Mode::Native: return sample_envmap(sample);
			case Sample::Mode::Cosine: return sample_cosine(sample);
			case Sample::Mode::Sphere: return sample_sphere(sample);
		}

		SLog(
			ELogLevel::EError,
			"This part of the code should never be executed. If it did, you likely passed an invalid Sample::Mode to Sample sample(Sample::Mode, Point2f&)."
		);
		return { };
	}

	void precompute()
	{
		Float result = 0.0;
		for (int y = 0; y < this->bitmap->getHeight(); ++y)
		{
			Float result_row = 0.0;
			for (int x = 0; x < this->bitmap->getWidth(); ++x)
			{
				Point2i pt(x, y);
				result_row += this->bitmap->getPixel(pt).getLuminance();
			}
			result += result_row;
			this->row_avgs.push_back(result_row / bitmap->getWidth());
		}
		SAssert(result != 0);
		this->bitmap_integral = result / (bitmap->getHeight() * bitmap->getWidth());
		this->precomputed = true;
	}

	/// Generates an empty bitmap with the same params as this one
	ref<Bitmap> gen_empty_bitmap()
	{
		Bitmap::EPixelFormat px_format = this->bitmap->getPixelFormat();
		Bitmap::EComponentFormat cmp_format = this->bitmap->getComponentFormat();
		Vector2i size = this->bitmap->getSize();
		std::size_t channels = this->bitmap->getChannelCount();

		return new Bitmap(px_format, cmp_format, size, channels, nullptr);
	}

	/// Noisifies the underlying bitmap via a custom Gaussian noise implementation
	void noisify(float noise_perc = 0.2f)
	{
		if (this->precomputed)
		{
			SLog(
				ELogLevel::EError,
				"An attempt was made to noisify after envmap precomputation. Make sure to noisify first and then precompute."
			);
		}

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
private:
	Sample sample_helper(Vector& dir)
	{
		/* Transform to (hemi)spherical coordinates */
		Float theta = std::acos(dir.z);
		Float phi = std::atan2(dir.y, dir.x);

		Point2f uv_norm(
			0.5f - phi * INV_TWOPI, // normalize phi into range [0.0, 1.0)
			theta * INV_PI 			// normalize theta into range [0.0, 0.5)
		);

		Point2i uv(
			uv_norm.x * this->bitmap->getWidth(),
			uv_norm.y * this->bitmap->getHeight()
		);

		/* UNCOMMENT FOR SAMPLE VISUALIZATION */
		sampled_points.push_back(std::pair<Point2i, Spectrum>(uv, this->bitmap->getPixel(uv)));

		Sample sample_data = {
			.value = this->bitmap->getPixel(uv).getLuminance(),
			.phi = phi + M_PI, // we want to move phi from [-pi, pi] to [0, 2pi]
			.theta = theta
		};

		/* Return found texel */
		return sample_data;
	}

	Sample sample_cosine(Point2f& sample)
	{
		Vector dir = warp::squareToCosineHemisphere(sample);
		return sample_helper(dir);
	}

	Sample sample_sphere(Point2f& sample)
	{
		Vector dir = warp::squareToUniformSphere(sample);
		return sample_helper(dir);
	}

	Sample sample_envmap(Point2f& sample)
	{
		Float sum_y = 0.0f;
		/* Iterate over rows until our sample is bigger than the respective avg. density */
		int y = 0;
		for (y = 0; y < this->row_avgs.size(); ++y)
		{
			sum_y += this->row_avgs.at(y) / this->bitmap_integral;
			if ((sum_y / this->bitmap->getHeight()) >= sample.y) break;
		}

		Float sum_x = 0.0f;
		/* Iterate over entries in row until our sample is bigger than the respective value */
		int x = 0;
		for (x = 0; x < this->bitmap->getWidth(); ++x)
		{
			Point2i pt(x, y);
			sum_x += this->bitmap->getPixel(pt).getLuminance() / this->row_avgs.at(y);
			if ((sum_x / this->bitmap->getWidth()) >= sample.x) break;
		}

		Point2i uv(x, y);

		/* UNCOMMENT FOR SAMPLE VISUALIZATION */
		sampled_points.push_back(std::pair<Point2i, Spectrum>(uv, this->bitmap->getPixel(uv)));
 
		Point2f uv_norm(
			(Float) uv.x / this->bitmap->getWidth(),
			(Float) uv.y / this->bitmap->getHeight()
		);

		Sample sample_data = {
			.value = this->bitmap->getPixel(uv).getLuminance(),
			.phi = 2 * M_PI * uv_norm.x,
			.theta = M_PI * uv_norm.y
		};
		
		return sample_data;
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