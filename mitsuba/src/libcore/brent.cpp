/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#include <mitsuba/core/brent.h>

/**
 * \file brent.cpp
 * \brief Brent's method nonlinear zero finder
 *
 * The implementation is transcribed from the Apache Commons
 * Java implementation.
 */
MTS_NAMESPACE_BEGIN

BrentSolver::Result BrentSolver::solve(const boost::function<float (float)> &f,
		float min, float max) const {
	// return the first endpoint if it is good enough
	float yMin = f(min);
	if (std::abs(yMin) <= m_absAccuracy)
		return Result(true, 0, min, yMin);

	// return the second endpoint if it is good enough
	float yMax = f(max);
	if (std::abs(yMax) <= m_absAccuracy)
		return Result(true, 0, max, yMax);

	float sign = yMin * yMax;
	if (sign > 0) {
		SLog(EWarn, "BrentSolver: Function values at the endpoints do not have different signs -- "
			"endpoints: [%f, %f], values: [%f, %f]", min, max, yMin, yMax);
		return Result(false, 0, 0, 0);
	} else {
		// solve using only the first endpoint as initial guess
		return solve(f, min, yMin, max, yMax, min, yMin);
	}
}

BrentSolver::Result BrentSolver::solve(const boost::function<float (float)> &f,
		float min, float max, float initial) const {
	if (initial < min || initial > max) {
		SLog(EWarn, "BrentSolver: Invalid interval: lower=%f, initial=%f, upper=%f",
			min, max, initial);
		return Result(false, 0, 0, 0);
	}

	// return the initial guess if it is good enough
	float yInitial = f(initial);
	if (std::abs(yInitial) <= m_absAccuracy)
		return Result(true, 0, initial, yInitial);

	// return the first endpoint if it is good enough
	float yMin = f(min);
	if (std::abs(yMin) <= m_absAccuracy)
		return Result(true, 0, min, yMin);

	// reduce interval if min and initial bracket the root
	if (yInitial * yMin < 0)
		return solve(f, min, yMin, initial, yInitial, min, yMin);

	// return the second endpoint if it is good enough
	float yMax = f(max);
	if (std::abs(yMax) <= m_absAccuracy)
		return Result(true, 0, max, yMax);

	// reduce interval if initial and max bracket the root
	if (yInitial * yMax < 0)
		return solve(f, initial, yInitial, max, yMax, initial, yInitial);

	SLog(EWarn, "BrentSolver: Function values at the endpoints do not have different signs -- "
		"endpoints: [%f, %f], values: [%f, %f]", min, max, yMin, yMax);

	return Result(false, 0, 0, 0);
}

BrentSolver::Result BrentSolver::solve(const boost::function<float (float)> &f,
			 float x0, float y0,
			 float x1, float y1,
			 float x2, float y2) const {
	float delta = x1 - x0;
	float oldDelta = delta;

	size_t i = 0;
	while (i < m_maxIterations) {
		if (std::abs(y2) < std::abs(y1)) {
			// use the bracket point if is better than last approximation
			x0 = x1;
			x1 = x2;
			x2 = x0;
			y0 = y1;
			y1 = y2;
			y2 = y0;
		}
		if (std::abs(y1) <= m_absAccuracy) {
			// Avoid division by very small values. Assume
			// the iteration has converged (the problem may
			// still be ill conditioned)
			return Result(true, i, x1, y1);
		}
		float dx = x2 - x1;
		float tolerance =
			std::max(m_relAccuracyPos * std::abs(x1), m_absAccuracyPos);

		if (std::abs(dx) <= tolerance)
			return Result(true, i, x1, y1);
		if ((std::abs(oldDelta) < tolerance) ||
				(std::abs(y0) <= std::abs(y1))) {
			// Force bisection.
			delta = (float) 0.5f * dx;
			oldDelta = delta;
		} else {
			float r3 = y1 / y0;
			float p;
			float p1;
			// the equality test (x0 == x2) is intentional,
			// it is part of the original Brent's method,
			// it should NOT be replaced by proximity test
			if (x0 == x2) {
				// Linear interpolation.
				p = dx * r3;
				p1 = 1 - r3;
			} else {
				// Inverse quadratic interpolation.
				float r1 = y0 / y2;
				float r2 = y1 / y2;
				p = r3 * (dx * r1 * (r1 - r2) - (x1 - x0) * (r2 - 1));
				p1 = (r1 - 1) * (r2 - 1) * (r3 - 1);
			}
			if (p > 0) {
				p1 = -p1;
			} else {
				p = -p;
			}
			if (2 * p >= (float) 1.5f * dx * p1 - std::abs(tolerance * p1) ||
					p >= std::abs((float) 0.5f * oldDelta * p1)) {
				// Inverse quadratic interpolation gives a value
				// in the wrong direction, or progress is slow.
				// Fall back to bisection.
				delta = (float) 0.5f * dx;
				oldDelta = delta;
			} else {
				oldDelta = delta;
				delta = p / p1;
			}
		}
		// Save old X1, Y1
		x0 = x1;
		y0 = y1;
		// Compute new X1, Y1
		if (std::abs(delta) > tolerance) {
			x1 = x1 + delta;
		} else if (dx > 0) {
			x1 = x1 + (float) 0.5f * tolerance;
		} else if (dx <= 0) {
			x1 = x1 - (float) 0.5f * tolerance;
		}
		y1 = f(x1);
		if ((y1 > 0) == (y2 > 0)) {
			x2 = x0;
			y2 = y0;
			delta = x1 - x0;
			oldDelta = delta;
		}
		i++;
	}
	SLog(EWarn, "BrentSolver: Maximum number of iterations (" SIZE_T_FMT ") exceeded!",
		m_maxIterations);
	return Result(false, i, x1, y1);
}

MTS_NAMESPACE_END
