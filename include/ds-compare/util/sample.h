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

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_UTIL_SAMPLE_H_ */