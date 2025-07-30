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
	struct TransformationPair {
		std::function<T(T)> f_org;
		std::function<T(T)> f_inv;

		TransformationPair(std::function<T(T)> org, std::function<T(T)> inv) : f_org(org), f_inv(inv) { };
	};

	template <typename T>
	static std::function<T(T)> f(Transformation t, bool inverse)
	{
		static TransformationPair<T> planar(
			[](T value) { return value; },
			[](T value) { return value; }
		);

		static TransformationPair<T> spherical(
			[](T value) {
				auto j = std::floor(value);
				auto pos = value - j;
				auto i = (2.0 * pos) - 1.0;
				return (std::acos(-i) * INV_PI) + j;
			},
			[](T value) {
				auto j = std::floor(value);
				auto pos = value - j;
				return (1.0 - std::cos(pos * M_PI)) * 0.5 + j;
			}
		);

		static TransformationPair<T> cosine(
			[](T value) { return value; },
			[](T value) { return value; }
		);

		// WIP. Not verified, yields weird results.
		/*static auto cosine = [](T value) {
			float i = std::fmod(value, 0.5 * M_PI);
			if (i < 0) i += 0.5 * M_PI;
			int j = std::floor(value / (0.5 * M_PI));
			return (j + 0.5 * (1 - std::pow(-1, j) * std::cos(2 * i)));
		};*/

		if (t == Planar)	return inverse ? planar.f_inv    : planar.f_org;
		if (t == Spherical) return inverse ? spherical.f_inv : spherical.f_org;
		if (t == Cosine)	return inverse ? spherical.f_inv : spherical.f_org;

		return nullptr;
	}
};

#endif /* __DSCOMPARE_UTIL_DSPARAMS_H_ */