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

#pragma once
#if !defined(__MITSUBA_CORE_SHVECTOR_H_)
#define __MITSUBA_CORE_SHVECTOR_H_

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/quad.h>
#include <Eigen/Core>

MTS_NAMESPACE_BEGIN

/* Precompute normalization coefficients for the first 10 bands */
#define SH_NORMTBL_SIZE 10

struct SHVector;

/**
 * \brief Stores the diagonal blocks of a spherical harmonic
 * rotation matrix
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct MTS_EXPORT_CORE SHRotation {
	typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Matrix;

	std::vector<Matrix> blocks;

	/// Construct a new rotation storage for the given number of bands
	inline SHRotation(int bands) : blocks(bands) {
		for (int i=0; i<bands; ++i) {
			int dim = 2*i+1;
			blocks[i] = Matrix(dim, dim);
		}
	}

	/**
	 * \brief Transform a coefficient vector and store the result into
	 * the given target vector.
	 *
	 * The source and target must have the same number of bands.
	 */
	void operator()(const SHVector &source, SHVector &target) const;
};

/**
 * \brief Stores a truncated real spherical harmonics representation of
 * an L2-integrable function.
 *
 * Also provides some other useful functionality, such as evaluation,
 * projection and rotation.
 *
 * The Mathematica equivalent of the basis functions implemented here is:
 *
 * \code
 * SphericalHarmonicQ[l_, m_, \[Theta]_, \[Phi]_] :=
 *   Piecewise[{
 *      {SphericalHarmonicY[l, m, \[Theta], \[Phi]], m == 0},
 *      {Sqrt[2]*Re[SphericalHarmonicY[l, m, \[Theta], \[Phi]]], m > 0},
 *      {Sqrt[2]*Im[SphericalHarmonicY[l, -m, \[Theta], \[Phi]]], m < 0}
 *  }]
 * \endcode
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct MTS_EXPORT_CORE SHVector {
public:
	/// Construct an invalid SH vector
	inline SHVector()
		: m_bands(0) {
	}

	/// Construct a new SH vector (initialized to zero)
	inline SHVector(int bands)
		: m_bands(bands), m_coeffs(bands*bands) {
		clear();
	}

	/// Unserialize a SH vector to a binary data stream
	SHVector(Stream *stream);

	/// Copy constructor
	inline SHVector(const SHVector &v) : m_bands(v.m_bands),
		m_coeffs(v.m_coeffs) {
	}

	/// Return the number of stored SH coefficient bands
	inline int getBands() const {
		return m_bands;
	}

	/// Serialize a SH vector to a binary data stream
	void serialize(Stream *stream) const;

	/// Get the energy per band
	inline float energy(int band) const {
		float result = 0;
		for (int m=-band; m<=band; ++m)
			result += std::abs(operator()(band,m));
		return result;
	}

	/// Assignment
	inline SHVector &operator=(const SHVector &v) {
		m_bands = v.m_bands;
		m_coeffs = v.m_coeffs;
		return *this;
	}

	/// Set all coefficients to zero
	inline void clear() {
		m_coeffs.setZero();
	}

	/// Component-wise addition
	inline SHVector& operator+=(const SHVector &v) {
		ptrdiff_t extendBy = v.m_coeffs.size() - m_coeffs.size();
		if (extendBy > 0) {
			m_coeffs.conservativeResize(v.m_coeffs.rows());
			m_coeffs.tail(extendBy).setZero();
			m_bands = v.m_bands;
		}
		m_coeffs.head(m_coeffs.size()) += v.m_coeffs.head(m_coeffs.size());
		return *this;
	}

	/// Component-wise addition
	inline SHVector operator+(const SHVector &v) const {
		SHVector vec(std::max(m_bands, v.m_bands));
		if (m_bands > v.m_bands) {
			vec.m_coeffs = m_coeffs;
			vec.m_coeffs.head(v.m_coeffs.size()) += v.m_coeffs;
		} else {
			vec.m_coeffs = v.m_coeffs;
			vec.m_coeffs.head(m_coeffs.size()) += m_coeffs;
		}
		return vec;
	}

	/// Component-wise subtraction
	inline SHVector& operator-=(const SHVector &v) {
		ptrdiff_t extendBy = v.m_coeffs.size() - m_coeffs.size();
		if (extendBy > 0) {
			m_coeffs.conservativeResize(v.m_coeffs.rows());
			m_coeffs.tail(extendBy).setZero();
			m_bands = v.m_bands;
		}
		m_coeffs.head(m_coeffs.size()) -= v.m_coeffs.head(m_coeffs.size());
		return *this;
	}

	/// Component-wise subtraction
	inline SHVector operator-(const SHVector &v) const {
		SHVector vec(std::max(m_bands, v.m_bands));
		if (m_bands > v.m_bands) {
			vec.m_coeffs = m_coeffs;
			vec.m_coeffs.head(v.m_coeffs.size()) -= v.m_coeffs;
		} else {
			vec.m_coeffs = -v.m_coeffs;
			vec.m_coeffs.head(m_coeffs.size()) += m_coeffs;
		}
		return vec;
	}

	/// Add a scalar multiple of another vector
	inline SHVector& madd(float f, const SHVector &v) {
		ptrdiff_t extendBy = v.m_coeffs.size() - m_coeffs.size();
		if (extendBy > 0) {
			m_coeffs.conservativeResize(v.m_coeffs.rows());
			m_coeffs.tail(extendBy).setZero();
			m_bands = v.m_bands;
		}
		m_coeffs.head(m_coeffs.size()) += v.m_coeffs.head(m_coeffs.size()) * f;

		return *this;
	}

	/// Scalar multiplication
	inline SHVector &operator*=(float f) {
		m_coeffs *= f;
		return *this;
	}

	/// Scalar multiplication
	inline SHVector operator*(float f) const {
		SHVector vec(m_bands);
		vec.m_coeffs = m_coeffs * f;
		return vec;
	}

	/// Scalar division
	inline SHVector &operator/=(float f) {
		m_coeffs *= (float) 1 / f;
		return *this;
	}

	/// Scalar division
	inline SHVector operator/(float f) const {
		SHVector vec(m_bands);
		vec.m_coeffs = m_coeffs * (1/f);
		return vec;
	}

	/// Negation operator
	inline SHVector operator-() const {
		SHVector vec(m_bands);
		vec.m_coeffs = -m_coeffs;
		return vec;
	}

	/// Access coefficient m (in {-l, ..., l}) on band l
	inline float &operator()(int l, int m) {
		return m_coeffs[l*(l+1) + m];
	}

	/// Access coefficient m (in {-l, ..., l}) on band l
	inline const float &operator()(int l, int m) const {
		return m_coeffs[l*(l+1) + m];
	}

	/// Evaluate for a direction given in spherical coordinates
	float eval(float theta, float phi) const;

	/// Evaluate for a direction given in Cartesian coordinates
	float eval(const Vector &v) const;

	/**
	 * \brief Evaluate for a direction given in spherical coordinates.
	 *
	 * This function is much faster but only works for azimuthally
	 * invariant functions
	 */
	float evalAzimuthallyInvariant(float theta, float phi) const;

	/**
	 * \brief Evaluate for a direction given in cartesian coordinates.
	 *
	 * This function is much faster but only works for azimuthally
	 * invariant functions
	 */
	float evalAzimuthallyInvariant(const Vector &v) const;

	/// Check if this function is azumuthally invariant
	bool isAzimuthallyInvariant() const;

	/// Equality comparison operator
	inline bool operator==(const SHVector &v) const {
		return m_bands == v.m_bands && m_coeffs == v.m_coeffs;
	}

	/// Equality comparison operator
	inline bool operator!=(const SHVector &v) const {
		return !operator==(v);
	}

	/// Dot product
	inline friend float dot(const SHVector &v1, const SHVector &v2);

	/// Normalize so that the represented function becomes a valid distribution
	void normalize();

	/// Compute the second spherical moment (analytic)
	Matrix3x3 mu2() const;

	/// Brute-force search for the minimum value over the sphere
	float findMinimum(int res) const;

	/// Add a constant value
	void addOffset(float value);

	/**
	 * \brief Convolve the SH representation with the supplied kernel.
	 *
	 * Based on the Funk-Hecke theorem -- the kernel must be rotationally
	 * symmetric around the Z-axis.
	 */
	void convolve(const SHVector &kernel);

	/// Project the given function onto a SH basis (using a 2D composite Simpson's rule)
	template<typename Functor> void project(const Functor &f, int res = 32) {
		SAssert(res % 2 == 0);
		/* Nested composite Simpson's rule */
		float hExt = M_PI / res,
		      hInt = (2*M_PI)/(res*2);

		for (int l=0; l<m_bands; ++l)
			for (int m=-l; m<=l; ++m)
				operator()(l,m) = 0;

		float *sinPhi = (float *) alloca(sizeof(float)*m_bands),
			  *cosPhi = (float *) alloca(sizeof(float)*m_bands);

		for (int i=0; i<=res; ++i) {
			float theta = hExt*i, cosTheta = std::cos(theta);
			float weightExt = (i & 1) ? 4.0f : 2.0f;
			if (i == 0 || i == res)
				weightExt = 1.0f;

			for (int j=0; j<=res*2; ++j) {
				float phi = hInt*j;
				float weightInt = (j & 1) ? 4.0f : 2.0f;
				if (j == 0 || j == 2*res)
					weightInt = 1.0f;

				for (int m=0; m<m_bands; ++m) {
					sinPhi[m] = std::sin((m+1)*phi);
					cosPhi[m] = std::cos((m+1)*phi);
				}

				float value = f(sphericalDirection(theta, phi))*std::sin(theta)
					* weightExt*weightInt;

				for (int l=0; l<m_bands; ++l) {
					for (int m=1; m<=l; ++m) {
						float L = legendreP(l, m, cosTheta) * normalization(l, m);
						operator()(l, -m) += value * SQRT_TWO * sinPhi[m-1] * L;
						operator()(l, m) += value * SQRT_TWO * cosPhi[m-1] * L;
					}

					operator()(l, 0) += value * legendreP(l, 0, cosTheta) * normalization(l, 0);
				}
			}
		}

		for (int l=0; l<m_bands; ++l)
			for (int m=-l; m<=l; ++m)
				operator()(l,m) *= hExt*hInt/9;
	}

	/// Compute the relative L2 error
	template<typename Functor> float l2Error(const Functor &f, int res = 32) const {
		SAssert(res % 2 == 0);
		/* Nested composite Simpson's rule */
		float hExt = M_PI / res,
		      hInt = (2*M_PI)/(res*2);
		float error = 0.0f, denom=0.0f;

		for (int i=0; i<=res; ++i) {
			float theta = hExt*i;
			float weightExt = (i & 1) ? 4.0f : 2.0f;
			if (i == 0 || i == res)
				weightExt = 1.0f;

			for (int j=0; j<=res*2; ++j) {
				float phi = hInt*j;
				float weightInt = (j & 1) ? 4.0f : 2.0f;
				if (j == 0 || j == 2*res)
					weightInt = 1.0f;

				float value1 = f(sphericalDirection(theta, phi));
				float value2 = eval(theta, phi);
				float diff = value1-value2;
				float weight = std::sin(theta)*weightInt*weightExt;

				error += diff*diff*weight;
				denom += value1*value1*weight;
			}
		}

		return error/denom;
	}

	/// Turn into a string representation
	std::string toString() const;

	/// Return a normalization coefficient
	inline static float normalization(int l, int m) {
		if (l < SH_NORMTBL_SIZE)
			return m_normalization[l*(l+1)/2 + m];
		else
			return computeNormalization(l, m);
	}

	/**
	 * \brief Recursively computes rotation matrices for each band of SH coefficients.
	 *
	 * Based on 'Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion'
	 * by Ivanic and Ruedenberg. The implemented tables follow the notation in
	 * 'Spherical Harmonic Lighting: The Gritty Details' by Robin Green.
	 */
	static void rotation(const Transform &t, SHRotation &rot);

	/// Precomputes normalization coefficients for the first few bands
	static void staticInitialization();

	/// Free the memory taken up by staticInitialization()
	static void staticShutdown();
protected:
	/// Helper function for rotation() -- computes a diagonal block based on the previous level
	static void rotationBlock(const SHRotation::Matrix &M1, const SHRotation::Matrix &Mp, SHRotation::Matrix &Mn);

	/// Compute a normalization coefficient
	static float computeNormalization(int l, int m);
private:
	int m_bands;
	Eigen::Matrix<float, Eigen::Dynamic, 1> m_coeffs;
	static float *m_normalization;
};

inline float dot(const SHVector &v1, const SHVector &v2) {
	const size_t size = std::min(v1.m_coeffs.size(), v2.m_coeffs.size());
	return v1.m_coeffs.head(size).dot(v2.m_coeffs.head(size));
}

/**
 * \brief Implementation of 'Importance Sampling Spherical Harmonics'
 * by W. Jarsz, N. Carr and H. W. Jensen (EUROGRAPHICS 2009)
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE SHSampler : public Object {
public:
	/**
	 * \brief Precompute a spherical harmonics sampler object
	 *
	 * \param bands Number of SH coefficient bands to support
	 * \param depth Number of recursive sample warping steps.
	 */
	SHSampler(int bands, int depth);

	/**
	 * \brief Warp a uniform sample in [0,1]^2 to one that is
	 * approximately proportional to the specified function.
	 *
	 * The resulting sample will have spherical coordinates
	 * [0,pi]x[0,2pi] and its actual PDF (which might be
	 * slightly different from the function evaluated at the
	 * sample, even if $f$ is a distribution) will be returned.
	 */
	float warp(const SHVector &f, Point2 &sample) const;

	/// Return information on the size of the precomputed tables
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~SHSampler();

	/* Index into the assoc. legendre polynomial table */
	inline int I(int l, int m) const { return l*(l+1)/2 + m; }

	/* Index into the phi table */
	inline int P(int m) const { return m + m_bands; }

	inline float lookupIntegral(int depth, int zBlock, int phiBlock, int l, int m) const {
		return -m_phiMap[depth][phiBlock][P(m)] * m_legendreMap[depth][zBlock][I(l, std::abs(m))];
	}

	/// Recursively compute assoc. legendre & phi integrals
	float *legendreIntegrals(float a, float b);
	float *phiIntegrals(float a, float b);

	/// Integrate a SH expansion over the specified mip-map region
	float integrate(int depth, int zBlock, int phiBlock, const SHVector &f) const;
protected:
	int m_bands;
	int m_depth;
	float ***m_phiMap;
	float ***m_legendreMap;
	int m_dataSize;
	float *m_normalization;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_SHVECTOR_H_ */
