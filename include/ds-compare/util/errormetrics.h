#pragma once

#if !defined(__DSCOMPARE_UTIL_ERRORMETRICS_H_)
#define __DSCOMPARE_UTIL_ERRORMETRICS_H_

#include <ds-compare/util/_definitions.h>
#include <ds-compare/util/envmap.h>

MTS_NAMESPACE_BEGIN

struct DS_COMPARE ErrorMetrics {
	static float MSE(const EnvironmentMap& em1, const EnvironmentMap& em2)
	{
		SAssert(em1.bitmap->getSize() == em2.bitmap->getSize());

		const Vector2i size = em1.bitmap->getSize();
		const int channels = em1.bitmap->getChannelCount();

		float err = 0.0f;
		for (int y = 0; y < size.y; ++y)
		{
			for (int x = 0; x < size.x; ++x)
			{
				const Point2i pt(x, y);
				Point3 em1_rgb = em1.get_pixel_rgb(pt);
				Point3 em2_rgb = em2.get_pixel_rgb(pt);

				// Calculate error at given pixel by squaring the error
				float sum = 0.0f;
				for (int i = 0; i < channels; ++i)
				{
					sum += (em1_rgb[i] - em2_rgb[i]) * (em1_rgb[i] - em2_rgb[i]);
				}

				err += sum;
			}
		}

		// Divide by the pixel count
		err /= em1.bitmap->getPixelCount() * channels;
		return err;
	}

	static float RMSE(const EnvironmentMap& em1, const EnvironmentMap& em2)
	{
		// Take the square root of MSE for RMSE
		return std::sqrt(ErrorMetrics::MSE(em1, em2));
	}

	static float MAE(const EnvironmentMap& em1, const EnvironmentMap& em2)
	{
		SAssert(em1.bitmap->getSize() == em2.bitmap->getSize());

		const Vector2i size = em1.bitmap->getSize();
		const int channels = em1.bitmap->getChannelCount();

		float err = 0.0f;
		for (int y = 0; y < size.y; ++y)
		{
			for (int x = 0; x < size.x; ++x)
			{
				const Point2i pt(x, y);
				Point3 em1_rgb = em1.get_pixel_rgb(pt);
				Point3 em2_rgb = em2.get_pixel_rgb(pt);

				// Calculate error at given pixel by taking the abs of the difference
				float sum = 0.0f;
				for (int i = 0; i < channels; ++i)
				{
					sum += std::abs(em1_rgb[i] - em2_rgb[i]);
				}

				err += sum;
			}
		}

		// Divide by the pixel count
		err /= em1.bitmap->getPixelCount() * channels;
		return err;
	}

	template <typename T, typename std::enable_if<std::is_arithmetic<T>::value>::type* = nullptr>
	static float MD(const std::vector<T>& reference, const Float gt_mean)
	{
		std::size_t sample_count = reference.size();
		SAssert(sample_count > 0);

		float sum = 0.0f;
		for (const auto& value : reference)
		{
			sum += std::abs(value - gt_mean);
		}

		// Uncomment for estimator mean vs. gt mean comparison in console
		/*Float mean = 0.0;
		for (const auto& observation : reference)
		{
			mean += observation;
		}
		mean /= reference.size();

		std::cout << "[Estimator / GT]" << mean << " / " << gt_mean << std::endl;
		std::cout << "MD: " << (sum / sample_count) << std::endl;*/

		return (sum / sample_count);
	}
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_UTIL_ERRORMETRICS_H_ */