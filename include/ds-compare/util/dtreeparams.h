#pragma once

#if !defined(__DSCOMPARE_UTIL_DTREEPARAMS_H_)
#define __DSCOMPARE_UTIL_DTREEPARAMS_H_

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

#endif /* __DSCOMPARE_UTIL_SAMPLE_H_ */