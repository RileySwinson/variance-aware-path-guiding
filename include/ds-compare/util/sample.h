#pragma once

#if !defined(__DSCOMPARE_UTIL_SAMPLE_H_)
#define __DSCOMPARE_UTIL_SAMPLE_H_

#include <ds-compare/util/_definitions.h>

MTS_NAMESPACE_BEGIN

struct DS_COMPARE Sample {
    Float value = 0;
	Float pdf = 0;

	Float theta = 0;
    Float phi = 0;

	enum Mode {
		Cosine,
		Native,
		Sphere
	};

	enum Strategy {
		Preprocess,
		Forward
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

	friend std::istream& operator>>(std::istream& in, Sample::Strategy& strategy)
	{
		std::string token;
		in >> token;

		if (token == "preprocess")
		{
			strategy = Sample::Strategy::Preprocess;
			return in;
		}
		if (token == "forward")
		{
			strategy = Sample::Strategy::Forward;
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

	std::string to_string()
	{
		return "value » " + std::to_string(this->value) + "\n"
			+ "pdf » " + std::to_string(this->pdf) + "\n"
			+ "φ » " + std::to_string(this->phi) + " | θ » " + std::to_string(this->theta);
	}

	bool is_valid()
	{
		if (!this->pdf || !std::isfinite(this->theta) || !std::isfinite(this->phi)) return false;

		bool phi_valid = (0 <= this->phi) && (this->phi < 2 * M_PI);
		bool theta_valid = (0 <= this->theta) && (this->theta < M_PI);

		return phi_valid && theta_valid;
	}
};

struct SampleStorage {
	SampleStorage(const Vector2i& map_dims) : map_dims(map_dims) { }

	void store(const Point2i& im_coords, Sample sample)
	{
		int pos = (im_coords.y * this->map_dims.x) + im_coords.x;
		this->umap_samples[pos].push_back(sample);
	}

	std::vector<Sample> obtain(const Point2i& im_coords)
	{
		int pos = (im_coords.y * this->map_dims.x) + im_coords.x;
		
		auto it = this->umap_samples.find(pos);
		if (it == this->umap_samples.end())
		{
			return std::vector<Sample>();
		}

		return it->second;
	}

	inline std::size_t max() const
	{
		return std::max_element(
			this->umap_samples.begin(), 
			this->umap_samples.end(),
			[](const auto& a, const auto& b) { 
				return a.second.size() < b.second.size();
			}
		)->second.size();
	}

	std::vector<Sample> to_flat(std::size_t init_size = 0) const
	{
		std::vector<Sample> flat_vector;
		flat_vector.reserve(init_size);

		for (const auto& observation : this->umap_samples)
		{
			flat_vector.insert(flat_vector.end(), observation.second.begin(), observation.second.end());
		}

		return flat_vector;
	}

	void write(const std::string& path)
	{
		std::ofstream output;
		output.open(path, std::ios::out);

		for (int y = 0; y < this->map_dims.y; ++y)
		{
			std::string data = "";

			for (int x = 0; x < this->map_dims.x; ++x)
			{
				Point2i pos(x, y);
				auto sample_count = obtain(pos).size();
				data += std::to_string(sample_count) + ",";
			}

			output << data;
			output << "\n";
		}

		output.close();
	}

private:
	std::unordered_map<uint32_t, std::vector<Sample>> umap_samples;
	const Vector2i map_dims;
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_UTIL_SAMPLE_H_ */