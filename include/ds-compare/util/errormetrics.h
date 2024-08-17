#pragma once

#if !defined(__DSCOMPARE_UTIL_ERRORMETRICS_H_)
#define __DSCOMPARE_UTIL_ERRORMETRICS_H_

#include <ds-compare/util/_definitions.h>

MTS_NAMESPACE_BEGIN

struct DS_COMPARE ErrorMetrics {
	static float MSE(const Bitmap& bm1, const Bitmap& bm2)
	{
		SAssert(bm1.getSize() == bm2.getSize());

		const int CHANNELS = 3;

		float err = 0.0f;
		for (int y = 0; y < bm1.getHeight(); ++y)
		{
			for (int x = 0; x < bm1.getWidth(); ++x)
			{
				Point2i pt(x, y);
				float bm1_rgb[CHANNELS];
				float bm2_rgb[CHANNELS];
				bm1.getPixel(pt).toLinearRGB(bm1_rgb[0], bm1_rgb[1], bm1_rgb[2]);
				bm2.getPixel(pt).toLinearRGB(bm2_rgb[0], bm2_rgb[1], bm2_rgb[2]);

				// Calculate error at given pixel by squaring the error
				float sum = 0.0f;
				for (int i = 0; i < CHANNELS; ++i)
				{
					sum += (bm1_rgb[i] - bm2_rgb[i]) * (bm1_rgb[i] - bm2_rgb[i]);
				}

				err += sum;
			}
		}

		// Divide by the pixel count
		err /= bm1.getPixelCount() * CHANNELS;
		return err;
	}

	static float RMSE(const Bitmap& bm1, const Bitmap& bm2)
	{
		// Take the square root of MSE for RMSE
		return std::sqrt(ErrorMetrics::MSE(bm1, bm2));
	}

	static float MAE(const Bitmap& bm1, const Bitmap& bm2)
	{
		SAssert(bm1.getSize() == bm2.getSize());

		const int CHANNELS = 3;

		float err = 0.0f;
		for (int y = 0; y < bm1.getHeight(); ++y)
		{
			for (int x = 0; x < bm1.getWidth(); ++x)
			{
				Point2i pt(x, y);
				float bm1_rgb[CHANNELS];
				float bm2_rgb[CHANNELS];
				bm1.getPixel(pt).toLinearRGB(bm1_rgb[0], bm1_rgb[1], bm1_rgb[2]);
				bm2.getPixel(pt).toLinearRGB(bm2_rgb[0], bm2_rgb[1], bm2_rgb[2]);

				// Calculate error at given pixel by taking the abs of the difference
				float sum = 0.0f;
				for (int i = 0; i < CHANNELS; ++i)
				{
					sum += std::abs(bm1_rgb[i] - bm2_rgb[i]);
				}

				err += sum;
			}
		}

		// Divide by the pixel count
		err /= bm1.getPixelCount() * CHANNELS;
		return err;
	}
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_UTIL_ERRORMETRICS_H_ */