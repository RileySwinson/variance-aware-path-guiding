#pragma once

#if !defined(__DSCOMPARE_UTIL_CONVERTER_H_)
#define __DSCOMPARE_UTIL_CONVERTER_H_

#include <ds-compare/util/_definitions.h>

MTS_NAMESPACE_BEGIN

struct DS_COMPARE Converter {
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

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_UTIL_CONVERTER_H_ */