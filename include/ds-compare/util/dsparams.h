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
				auto i = std::floor(value);
				auto f = value - i;
				auto pos = (2.0 * f) - 1.0;
				return (std::acos(-pos) * INV_PI) + i;
			},
			[](T value) {
				auto i = std::floor(value);
				auto f = value - i;
				return (1.0 - std::cos(f * M_PI)) * 0.5 + i;
			}
		);

		static TransformationPair<T> cosine(
			[](T value) {
				auto i = std::floor(value);
				auto sgn = (i % 2 == 0) ? 1 : -1;
				return 0.5 * (1 - sgn * std::cos(value * M_PI)) + i;
			},
			[](T value) {
				auto i = std::floor(value);
				auto f = value - i;
				auto pos = 1.0 - (2.0 * f);
				return (std::acos(pos) / M_PI) + i;
			}
		);

		if (t == Planar)	return inverse ? planar.f_inv    : planar.f_org;
		if (t == Spherical) return inverse ? spherical.f_inv : spherical.f_org;
		if (t == Cosine)	return inverse ? cosine.f_inv    : cosine.f_org;

		return nullptr;
	}
};

#endif /* __DSCOMPARE_UTIL_DSPARAMS_H_ */