#pragma once

#if !defined(__DSCOMPARE_UTIL_ENVMAP_H_)
#define __DSCOMPARE_UTIL_ENVMAP_H_

#include <ds-compare/util/_definitions.h>

MTS_NAMESPACE_BEGIN

struct EnvironmentMap;
typedef boost::optional<EnvironmentMap> OptionalEnvMap;

struct DS_COMPARE EnvironmentMap {
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

		EnvironmentMap::validate(envmap);
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
				const Point2i pt(x, y);
				result_row += get_pixel_luminance(pt);
			}
			result += result_row;
			this->row_avgs.push_back(result_row / bitmap->getWidth());
		}

		SAssert(result != 0);
		
		this->bitmap_integral = result / bitmap->getPixelCount();
		this->precomputed = true;
	}

	/// Get rgb color information at the specified position (x, y).
	inline Point3 get_pixel_rgb(const Point2i& pos) const
	{
		const Vector2i& size = this->bitmap->getSize();
		const int channels = this->bitmap->getChannelCount();

		SAssertEx(pos.x >= 0 && pos.x < size.x && pos.y >= 0 && pos.y < size.y, "EnvironmentMap::get_pixel_rgb(): out of bounds!");
		SAssertEx(channels == 3, "Envmap should have 3 channels!");
		SAssertEx(this->bitmap->getBytesPerComponent() == 4, "Envmap isn't float32!");

		size_t offset = ((size_t) pos.x + size.x * (size_t) pos.y) * channels;
		const float* data = this->bitmap->getFloat32Data();

		return Point3(
			data[offset],
			data[offset + 1],
			data[offset + 2]
		);
	}
	
	/// Get the luminance at a specified position (x, y).
	inline Float get_pixel_luminance(const Point2i& pos) const
	{
		Point3 rgb = get_pixel_rgb(pos);
		return rgb[0] * 0.212671f + rgb[1] * 0.715160f + rgb[2] * 0.072169f;
	}

	/// Set rgb color information at the specified position (x, y).
	inline void set_pixel_rgb(const Point2i& pos, const Point3& rgb)
	{
		const Vector2i& size = this->bitmap->getSize();
		const int channels = this->bitmap->getChannelCount();

		SAssertEx(pos.x >= 0 && pos.x < size.x && pos.y >= 0 && pos.y < size.y, "EnvironmentMap::get_pixel_rgb(): out of bounds!");
		SAssertEx(channels == 3, "Envmap should have 3 channels!");
		SAssertEx(this->bitmap->getBytesPerComponent() == 4, "Envmap isn't float32!");

		size_t offset = ((size_t) pos.x + size.x * (size_t) pos.y) * channels;
		float* data = this->bitmap->getFloat32Data();

		data[offset]	 = (float) rgb[0];
		data[offset + 1] = (float) rgb[1];
		data[offset + 2] = (float) rgb[2];
	}

	/// Generates an envmap with the same params as this one
	EnvironmentMap deep_copy(bool empty = false)
	{
		const Bitmap::EPixelFormat px_format = this->bitmap->getPixelFormat();
		const Bitmap::EComponentFormat cmp_format = this->bitmap->getComponentFormat();
		const Vector2i size = this->bitmap->getSize();
		const std::size_t channels = this->bitmap->getChannelCount();

		EnvironmentMap envmap;
		envmap.bitmap = new Bitmap(px_format, cmp_format, size, channels, NULL);

		for (int y = 0; y < size.y; ++y)
		{
			for (int x = 0; x < size.x; ++x)
			{
				Point2i pt(x, y);

				Point3 px(0, 0, 0);
				if (!empty) px = get_pixel_rgb(pt);

				envmap.set_pixel_rgb(pt, px);
			}
		}
		
		return envmap;
	}

	EnvironmentMap& map(std::function<void(Point2i coords, Point3& px)> F)
	{
		Vector2i size = this->bitmap->getSize();
		for (int y = 0; y < size.y; ++y)
		{
			for (int x = 0; x < size.x; ++x)
			{
				const Point2i pt(x, y);
				Point3 px = get_pixel_rgb(pt);

				F(pt, px); // do something with px...

				set_pixel_rgb(pt, px);
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

			Point3 px = get_pixel_rgb(pt);
			px /= max;

			set_pixel_rgb(pt, px);
		}

		return *this;
	}

	/// Calculate the weighted mean (w.r.t. solid angle) of the entire envmap
	Float mean()
	{
		Vector2i size = this->bitmap->getSize();

		if (size.x == 0 || size.y == 0)
		{
			return ((Float) 0);
		}

		Float sum_luminance = 0;
		Float sum_weights = 0;
		for (int y = 0; y < size.y; ++y)
		{
			Point2i temp(0, y);
			Point2 uv = Converter::image_to_uv(temp, size);
			Float theta = Converter::uv_to_spherical(uv).y;

			Float weight = std::max((Float) 0, std::sin(theta));

			for (int x = 0; x < size.x; ++x)
			{
				Point2i px(x, y);
				sum_luminance += get_pixel_luminance(px) * weight;
				sum_weights += weight;
			}
		}

		return (sum_luminance / sum_weights);
	}

	/// Calculate the weighted mean based on specific samples
	Float mean(std::vector<Sample> samples)
	{
		Vector2i size = this->bitmap->getSize();

		if (size.x == 0 || size.y == 0)
		{
			return ((Float) 0);
		}

		Float sum_luminance = 0;
		Float sum_weights = 0;
		for (const auto& sample : samples)
		{
			Float weight = std::max((Float) 0, std::sin(sample.theta));

			Point2 spherical(sample.phi, sample.theta);
			Point2 uv = Converter::spherical_to_uv(spherical);
			Point2i px = Converter::uv_to_image(uv, size);

			sum_luminance += get_pixel_luminance(px) * weight;
			sum_weights += weight;
		}

		return (sum_luminance / sum_weights);
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

			Float px_lum = get_pixel_luminance(uv);
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

				const auto pt = Point2i(x, y);
				Point3 px = get_pixel_rgb(pt);
				
				float noise = random->nextFloat() * STD_DEV + MEAN;
				
				for (int channel = 0; channel < 3; ++channel)
				{
					px[channel] -= noise; 
					if (px[channel] < 0.0f) px[channel] = 0.0f;
				}

				set_pixel_rgb(pt, px);
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
			.value = get_pixel_luminance(uv),
			.pdf = pdf(mode, uv_norm),
			.theta = theta,
			.phi = phi
		};

		// Hacky, but we need to ensure that pdf != 0, which happens in infinitesimal cases in our envmaps like uv.u = 0.
		if (sample_data.pdf == 0) sample_data.pdf = Epsilon;

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
			const Point2i pt(x, y);
			sum_x += get_pixel_luminance(pt) / this->row_avgs.at(y);
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
			.value = get_pixel_luminance(uv),
			.pdf = pdf(Sample::Mode::Native, uv_norm),
			.theta = spherical.y,
			.phi = spherical.x
		};
		
		return sample_data;
	}

	static void validate(EnvironmentMap& envmap)
	{
		const Bitmap::EPixelFormat pf = envmap.bitmap->getPixelFormat();
		const Bitmap::EComponentFormat cf = envmap.bitmap->getComponentFormat();

		if (pf == Bitmap::EPixelFormat::ERGB && cf == Bitmap::EComponentFormat::EFloat32) return;

		const Bitmap::EPixelFormat px_format = Bitmap::EPixelFormat::ERGB;
		const Bitmap::EComponentFormat cmp_format = Bitmap::EComponentFormat::EFloat32;
		const Vector2i size = envmap.bitmap->getSize();
		const std::size_t channels = 3;

		ref<Bitmap> bm = new Bitmap(px_format, cmp_format, size, channels, NULL);
		envmap.bitmap->convert(bm);
		envmap.bitmap = bm;

		SAssert(envmap.bitmap->getChannelCount() == 3);
		SAssert(envmap.bitmap->getBytesPerComponent() == 4);
	}
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_UTIL_ENVMAP_H_ */