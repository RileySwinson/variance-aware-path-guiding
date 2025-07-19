#pragma once

#if !defined(__DSCOMPARE_UTIL_CONVERTER_H_)
#define __DSCOMPARE_UTIL_CONVERTER_H_

#include <ds-compare/util/_definitions.h>

MTS_NAMESPACE_BEGIN

struct DS_COMPARE Converter {
private:
	template<typename T>
	class Mapper
	{
	public:
		explicit Mapper(T value) : value(value) {};

		Mapper<T>& from(const std::pair<T, T>& from_range)
		{
			this->from_range = from_range;
			return *this;
		}

		T to(const std::pair<T, T>& to_range)
		{
			if (this->from_range.first == this->from_range.second)
			{
				return to_range.first;
			}

			T slope = (to_range.second - to_range.first) / (T)(this->from_range.second - this->from_range.first);
			return (slope * (this->value - this->from_range.first) + to_range.first);
		}

	private:
		T value;
		std::pair<T, T> from_range;
	};

public:
	/// Converts spherical coordinates [phi, theta] in the domain [0, 2pi) x [0, pi] to uv coordinates in the domain [0, 1)^2
	static Point2 spherical_to_uv(const Point2& spherical)
	{
		Float u = spherical.x * INV_TWOPI;
		Float v = spherical.y * INV_PI;

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

	/// Converts uv coordinates in the domain [0, 1)^2 to image coordinates (e.g., 1920x1080)
	static Point2i uv_to_image(const Point2& uv, const Vector2i& dims)
	{
		return Point2i(
			boost::algorithm::clamp(uv.x * dims.x, 0, dims.x - 1),
			boost::algorithm::clamp(uv.y * dims.y, 0, dims.y - 1)
		);
	}

	/// Converts image coordinates to uv coordinates in the domain [0, 1)^2
	static Point2 image_to_uv(const Point2i& coords, const Vector2i& dims)
	{
		Float u = coords.x / (Float) dims.x;
		Float v = coords.y / (Float) dims.y;

		return Point2(
			boost::algorithm::clamp(u, 0, 1 - Epsilon),
			boost::algorithm::clamp(v, 0, 1 - Epsilon)
		);
	}

	/// Maps a value linearly from one range to another range
	template<typename T>
	static Mapper<T> map(T value)
	{
		return Mapper<T>(value);
	}
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_UTIL_CONVERTER_H_ */