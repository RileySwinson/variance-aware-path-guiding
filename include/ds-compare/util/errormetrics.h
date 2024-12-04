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

	static float PSNR(const EnvironmentMap& em1, const EnvironmentMap& em2)
	{
		// We assume the image is normalized (for now)!
		const double MAX = 1.0;
		return 20 * std::log10(MAX) - 10 * std::log10(ErrorMetrics::MSE(em1, em2));
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

	/*static float SSIM(const Bitmap& bm1, const Bitmap& bm2)
	{
		SAssert(bm1.getSize() == bm2.getSize());

		const int CHANNELS = 3;

		// Constants
		const double k1 = 0.01;
		const double k2 = 0.03;
		const double L = std::pow(2, bm1.getBytesPerPixel() * 8) - 1;
		
		const double c1 = (k1 * L) * (k1 * L);
		const double c2 = (k2 * L) * (k2 * L);

		// Create gaussian kernel
		auto get_gaussian = [](int N, double sigma)
		{
			double r, s = 2.0 * sigma * sigma;
			double sum = 0.0;
			std::vector<double> kernel(N * N);

			int half_N = N / 2;
			for (int y = -half_N; y <= half_N; ++y)
			{
				for (int x = -half_N; x <= half_N; ++x)
				{
					r = std::sqrt(x * x + y * y);
					auto res = (std::exp(-(r * r) / s)) / (M_PI * s);

					int i = (y + half_N) * N + (x + half_N);
					kernel.at(i) = res;
					sum += res;
				}
			}

			for (double& v : kernel) v /= sum;

			return kernel;
		};

		const int N = 11;
		auto gaussian = get_gaussian(N, 1.5);

		// Convolve 
		double err = 0.0f;
		for (int y = 0; y < bm1.getHeight(); ++y)
		{
			for (int x = 0; x < bm1.getWidth(); ++x)
			{
				double bm1_mean[CHANNELS];
				double bm2_mean[CHANNELS];
				double bm1_var[CHANNELS];
				double bm2_var[CHANNELS];
				double covar[CHANNELS];

				// Calculate means and first part of variance
				for (int i = 0; i < gaussian.size(); ++i)
				{
					Point2i kernel_pos(i % N, i / N);

					Point2i image_pos(
						x + (kernel_pos.x - (N / 2)),
						y + (kernel_pos.y - (N / 2))
					);

					if (image_pos.x < 0 || image_pos.x >= bm1.getWidth()) // ... wrap around and sample on the other side of the image
					{
						image_pos.x = (bm1.getWidth() + image_pos.x) % bm1.getWidth();
					}

					if (image_pos.y < 0) // mirror
					{
						image_pos.y = std::abs(image_pos.y);
					}
					else if (image_pos.y >= bm1.getHeight())
					{
						image_pos.y = 2 * bm1.getHeight() - image_pos.y - 1;
					}

					Spectrum px1 = bm1.getPixel(image_pos);
					Spectrum px2 = bm2.getPixel(image_pos);

					for (char c = 0; c < CHANNELS; ++c)
					{
						double val1 = px1[c] * gaussian.at(i);
						double val2 = px2[c] * gaussian.at(i);

						bm1_mean[c] += val1;
						bm2_mean[c] += val2;

						bm1_var[c] += val1 * px1[c];
						bm2_var[c] += val2 * px2[c];

						covar[c] += val1 * val2;
					}
				}

				double local_ssim = 0.0;
				for (char c = 0; c < CHANNELS; ++c)
				{
					// Subtract by mean squared for variance and by product of means for covariance
					bm1_var[c] -= bm1_mean[c] * bm1_mean[c];
					bm2_var[c] -= bm2_mean[c] * bm2_mean[c];
					covar[c] -= bm1_mean[c] * bm2_mean[c];

					local_ssim += ((2 * bm1_mean[c] * bm2_mean[c]) * (2 * covar[c] + c2)) /
						(((bm1_mean[c] * bm1_mean[c]) + (bm2_mean[c] * bm2_mean[c]) + c1) * ((bm1_var[c] * bm1_var[c]) + (bm2_var[c] * bm2_var[c]) + c2));
				}

				err += local_ssim / CHANNELS;
			}
		}

		return err / (float) bm1.getPixelCount();
	}*/

	template <typename T, typename std::enable_if<std::is_arithmetic<T>::value>::type* = nullptr>
	static float MD(const std::vector<T>& reference, const std::vector<T>& observed)
	{
		SAssert(reference.size() == observed.size());

		float mean_deviation = 0.0f;

		int samples = reference.size();
		for (int i = 0; i < samples; ++i)
		{
			mean_deviation += std::abs(reference.at(i) - observed.at(i));
		}

		return (mean_deviation / samples);
	}
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_UTIL_ERRORMETRICS_H_ */