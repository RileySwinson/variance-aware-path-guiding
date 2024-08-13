#pragma once 

#if !defined(__MITSUBA_DS_UTIL_H_)
#define __MITSUBA_DS_UTIL_H_

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/random.h>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/clamp.hpp>

#include <array>
#include <functional>
#include <memory>
#include <map>
#include <random>

MTS_NAMESPACE_BEGIN

struct Sample;
struct EnvironmentMap;
struct ErrorMetrics;
struct DataStructure;

typedef boost::optional<EnvironmentMap> OptionalEnvMap;

struct MTS_EXPORT_CORE Sample {
    Float value = 0;
	Float pdf = 0;

	Float theta = 0;
    Float phi = 0;

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

	void noisify(Float mean = 1, Float stddev = 0.4)
	{
		static std::random_device rd;
		static std::mt19937 gen(rd());

		std::normal_distribution<Float> nd(mean, stddev);
		Float noise_factor = nd(gen);
		
		this->value = boost::algorithm::clamp(this->value * noise_factor, 0, this->value);
	}
};

struct MTS_EXPORT_CORE Converter {
	/// Converts spherical coordinates [phi, theta] in the domain [0, 2pi) x [0, pi] to uv coordinates in the domain [0, 1)^2
	static Point2 spherical_to_uv(const Point2& sphere)
	{
		Float u = sphere.x * INV_TWOPI;
		Float v = sphere.y * INV_PI;

		return Point2(
			boost::algorithm::clamp(u, 0, 1 - Epsilon),
			boost::algorithm::clamp(v, 0, 1 - Epsilon)
		);
	}

	/// Converts uv coordinates in the domain [0, 1)^2 to spherical coordinates [phi, theta]
	static Point2 uv_to_spherical(const Point2& uv)
	{
		Float phi = 2 * M_PI * uv.x;
		Float theta = M_PI * uv.y;

		return Point2(
			boost::algorithm::clamp(phi, 0, 2 * M_PI - Epsilon),
			boost::algorithm::clamp(theta, 0, M_PI - Epsilon)
		);
	}
};

struct MTS_EXPORT_CORE EnvironmentMap {
	ref<Bitmap> bitmap;
	std::array<std::string, 2> path;
	bool precomputed = false;

	Float bitmap_integral;
	std::vector<Float> row_avgs;

	static OptionalEnvMap fetch(const boost::filesystem::path& path)
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
			.path = { path.parent_path().filename().string(), path.stem().string() }
		};

		return envmap;
	}

	Sample sample(const Sample::Mode mode, const Point2& sample)
	{
		switch(mode)
		{
			case Sample::Mode::Native: return sample_envmap(sample);
			case Sample::Mode::Cosine: return sample_cosine(sample);
			case Sample::Mode::Sphere: return sample_sphere(sample);
		}

		SLog(
			ELogLevel::EError,
			"This part of the code should never be executed. If it did, you likely passed an invalid Sample::Mode to Sample sample(Sample::Mode, Point2&)."
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

		this->bitmap_integral = result / bitmap->getPixelCount();
		this->precomputed = true;
	}

	/// Generates an envmap with the same params as this one
	EnvironmentMap deep_copy(bool empty = false)
	{
		Bitmap::EPixelFormat px_format = this->bitmap->getPixelFormat();
		Bitmap::EComponentFormat cmp_format = this->bitmap->getComponentFormat();
		Vector2i size = this->bitmap->getSize();
		std::size_t channels = this->bitmap->getChannelCount();

		ref<Bitmap> bm = new Bitmap(px_format, cmp_format, size, channels, NULL);
		for (int y = 0; y < size.y; ++y)
		{
			for (int x = 0; x < size.x; ++x)
			{
				Point2i pt(x, y);
				Spectrum px = bm->getPixel(pt);

				for (int c = 0; c < channels; ++c)
					px[c] = empty ? 0 : this->bitmap->getPixel(pt)[c];

				bm->setPixel(pt, px);
			}
		}

		EnvironmentMap envmap;
		envmap.bitmap = bm;
		
		return envmap;
	}

	EnvironmentMap& map(std::function<void(Point2i coords, Spectrum& px)> F)
	{
		Vector2i size = this->bitmap->getSize();
		for (int y = 0; y < size.y; ++y)
		{
			for (int x = 0; x < size.x; ++x)
			{
				Point2i pt(x, y);
				Spectrum px = this->bitmap->getPixel(pt);

				F(pt, px); // do something with px...

				this->bitmap->setPixel(pt, px);
			}
		}

		return *this;
	}

	EnvironmentMap& normalize(Float max)
	{
		int y = 0;
		for (std::size_t v_i = 0; v_i < this->bitmap->getPixelCount(); ++v_i)
		{
			if (v_i != 0 && (v_i % this->bitmap->getWidth()) == 0) y += 1;

			Point2i pt(v_i % this->bitmap->getWidth(), y);
			Spectrum px = this->bitmap->getPixel(pt);

			for (int c = 0; c < this->bitmap->getChannelCount(); ++c)
				px[c] /= max;

			this->bitmap->setPixel(pt, px);
		}

		return *this;
	}

	void write(std::string path)
	{
		this->bitmap->write(Bitmap::EFileFormat::EOpenEXR, path);
	}

	/// Obtain the density value at a given position [x, y] in the domain [0, 1)^2
	Float pdf(Sample::Mode mode, Point2& pos) const
	{
		SAssert(pos.x >= 0 && pos.x < 1 && pos.y >= 0 && pos.y < 1);

		if (mode == Sample::Mode::Sphere)
		{
			Point2 sphere_coords = Converter::uv_to_spherical(pos);
			return INV_FOURPI * std::sin(sphere_coords.y);
		}

		if (mode == Sample::Mode::Cosine)
		{
			Point2 cosine_coords = Converter::uv_to_spherical(pos);
			return INV_PI * std::cos(cosine_coords.y) * std::sin(cosine_coords.y);
		}

		if (mode == Sample::Mode::Native)
		{
			Point2i uv(
				pos.x * this->bitmap->getWidth(),
				pos.y * this->bitmap->getHeight()
			);

			Float px_lum = this->bitmap->getPixel(uv).getLuminance();
			Float total_lum = this->bitmap_integral * this->bitmap->getPixelCount();
			return (px_lum / total_lum);
		}
		
		return 0;
	}

	/// Noisifies the underlying bitmap via a custom Gaussian noise implementation
	void noisify(const float noise_perc = 0.2f)
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
	Sample sample_helper(const Vector& dir, const Sample::Mode mode)
	{
		/* Transform to (hemi)spherical coordinates */
		Float theta = std::acos(dir.z);
		Float phi = std::atan2(dir.y, dir.x);
		if (phi < 0) phi += 2 * M_PI;

		Point2 uv_norm = Converter::spherical_to_uv(Point2(phi, theta));

		Point2i uv(
			uv_norm.x * this->bitmap->getWidth(),
			uv_norm.y * this->bitmap->getHeight()
		);

		Sample sample_data = {
			.value = this->bitmap->getPixel(uv).getLuminance(),
			.pdf = pdf(mode, uv_norm),
			.theta = theta,
			.phi = phi
		};

		/* Return found texel */
		return sample_data;
	}

	Sample sample_cosine(const Point2& sample)
	{
		Vector dir = warp::squareToCosineHemisphere(sample);
		return sample_helper(dir, Sample::Mode::Cosine);
	}

	Sample sample_sphere(const Point2& sample)
	{
		Vector dir = warp::squareToUniformSphere(sample);
		return sample_helper(dir, Sample::Mode::Sphere);
	}

	Sample sample_envmap(const Point2& sample)
	{
		Float sum_y = 0;
		/* Iterate over rows until our sample is bigger than the respective avg. density */
		int y = 0;
		for (y = 0; y < this->row_avgs.size(); ++y)
		{
			sum_y += this->row_avgs.at(y) / this->bitmap_integral;
			if ((sum_y / this->bitmap->getHeight()) >= sample.y) break;
		}

		// In case (sum_y / height) is smaller than sample.y (happens with samples *very* close to 1), we just subtract by 1.
		// As it would land on the very last pixel anyways, we're not introducing any bias here.
		if (y == this->bitmap->getHeight()) y -= 1;

		Float sum_x = 0;
		/* Iterate over entries in row until our sample is bigger than the respective value */
		int x = 0;
		for (x = 0; x < this->bitmap->getWidth(); ++x)
		{
			Point2i pt(x, y);
			sum_x += this->bitmap->getPixel(pt).getLuminance() / this->row_avgs.at(y);
			if ((sum_x / this->bitmap->getWidth()) >= sample.x) break;
		}

		if (x == this->bitmap->getWidth()) x -= 1;

		Point2i uv(x, y);
		Point2 uv_norm(
			(Float) uv.x / this->bitmap->getWidth(),
			(Float) uv.y / this->bitmap->getHeight()
		);

		Point2 spherical = Converter::uv_to_spherical(uv_norm);

		Sample sample_data = {
			.value = this->bitmap->getPixel(uv).getLuminance(),
			.pdf = pdf(Sample::Mode::Native, uv_norm),
			.theta = spherical.y,
			.phi = spherical.x
		};
		
		return sample_data;
	}
};

struct MTS_EXPORT_CORE ErrorMetrics {
	static float MSE(const Bitmap& bm1, const Bitmap& bm2)
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

	static float RMSE(const Bitmap& bm1, const Bitmap& bm2)
	{
		// Take the square root of MSE for RMSE
		return std::sqrt(ErrorMetrics::MSE(bm1, bm2));
	}

	static float MAE(const Bitmap& bm1, const Bitmap& bm2)
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

struct MTS_EXPORT_CORE StatTrak { /* TODO */ };

MTS_NAMESPACE_END

#endif /* __MITSUBA_DS_UTIL_H_ */