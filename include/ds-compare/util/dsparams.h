#pragma once

#if !defined(__DSCOMPARE_UTIL_DSPARAMS_H_)
#define __DSCOMPARE_UTIL_DSPARAMS_H_

#include <ds-compare/util/_definitions.h>

struct DS_COMPARE DTreeParams {
	enum EBsdfSamplingFractionLoss {
		ENone,
		EKL,
		EVariance,
	};

	enum EDirectionalFilter {
		ENearest,
		EBox,
	};

	friend std::istream& operator>>(std::istream& in, EBsdfSamplingFractionLoss& mode)
	{
		std::string token;
		in >> token;

		if (token == "none")
		{
			mode = EBsdfSamplingFractionLoss::ENone;
			return in;
		}
		if (token == "kl")
		{
			mode = EBsdfSamplingFractionLoss::EKL;
			return in;
		}
		if (token == "var")
		{
			mode = EBsdfSamplingFractionLoss::EVariance;
			return in;
		}

		in.setstate(std::ios_base::failbit);
		return in;
	};

	friend std::istream& operator>>(std::istream& in, EDirectionalFilter& mode)
	{
		std::string token;
		in >> token;

		if (token == "nearest")
		{
			mode = EDirectionalFilter::ENearest;
			return in;
		}
		if (token == "box")
		{
			mode = EDirectionalFilter::EBox;
			return in;
		}

		in.setstate(std::ios_base::failbit);
		return in;
	};
};

struct DS_COMPARE TCParams {
	/**
	 * @brief Used to specify how values are transformed to account for geometrical (i.e., spherical, cosine, ...) properties.
	 */
	enum Transformation {
		Planar,
		Spherical,
		Cosine,
		// Feel free to add more here. Make sure to adjust the function below accordingly.
	};

	friend std::istream& operator>>(std::istream& in, Transformation& mode)
	{
		std::string token;
		in >> token;

		if (token == "planar")
		{
			mode = Transformation::Planar;
			return in;
		}
		if (token == "spherical")
		{
			mode = Transformation::Spherical;
			return in;
		}
		if (token == "cosine")
		{
			mode = Transformation::Cosine;
			return in;
		}

		in.setstate(std::ios_base::failbit);
		return in;
	}

	template <typename T>
	static std::function<T(T)> f(Transformation t)
	{
		static auto planar = [](T value) {
			return value;
		};

		static auto spherical = [](T value) {
			auto i = std::fmod(value + 1, 2.0) - 1.0;
			auto j = std::floor((value + 1) * 0.5);
			return (j * M_PI) + std::asin(i);
		};

		static auto cosine = [](T value) {
			// WIP. Not verified, yields weird results.
			float i = std::fmod(value, 0.5 * M_PI);
			if (i < 0) i += 0.5 * M_PI;
			int j = std::floor(value / (0.5 * M_PI));
			return (j + 0.5 * (1 - std::pow(-1, j) * std::cos(2 * i)));
		};

		if (t == Planar)	return planar;
		if (t == Spherical) return spherical;
		if (t == Cosine)	return cosine;

		return nullptr;
	}
};

#endif /* __DSCOMPARE_UTIL_DSPARAMS_H_ */