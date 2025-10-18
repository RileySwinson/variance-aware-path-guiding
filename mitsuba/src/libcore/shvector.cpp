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

#include <mitsuba/core/shvector.h>
#include <mitsuba/core/transform.h>
#include <boost/math/special_functions/factorials.hpp>

MTS_NAMESPACE_BEGIN

float *SHVector::m_normalization = NULL;

SHVector::SHVector(Stream *stream) {
	m_bands = stream->readInt();
	unsigned int size = m_bands*m_bands;
	m_coeffs.resize(size);
	for (size_t i=0; i<size; ++i)
		m_coeffs[i] = stream->readfloat();
}

void SHVector::serialize(Stream *stream) const {
	stream->writeInt(m_bands);
	for (size_t i=0; i<(size_t) m_coeffs.size(); ++i)
		stream->writefloat(m_coeffs[i]);
}

bool SHVector::isAzimuthallyInvariant() const {
	for (int l=0; l<m_bands; ++l) {
		for (int m=1; m<=l; ++m) {
			if (std::abs(operator()(l, -m)) > Epsilon
			 || std::abs(operator()(l, m)) > Epsilon)
				return false;
		}
	}
	return true;
}

float SHVector::eval(float theta, float phi) const {
	float result = 0;
	float cosTheta = std::cos(theta);
	float *sinPhi = (float *) alloca(sizeof(float)*m_bands),
		  *cosPhi = (float *) alloca(sizeof(float)*m_bands);

	for (int m=0; m<m_bands; ++m) {
		sinPhi[m] = std::sin((m+1) * phi);
		cosPhi[m] = std::cos((m+1) * phi);
	}

	for (int l=0; l<m_bands; ++l) {
		for (int m=1; m<=l; ++m) {
			float L = legendreP(l, m, cosTheta) * normalization(l, m);
			result += operator()(l, -m) * SQRT_TWO * sinPhi[m-1] * L;
			result += operator()(l, m)  * SQRT_TWO * cosPhi[m-1] * L;
		}

		result += operator()(l, 0) * legendreP(l, 0, cosTheta) * normalization(l, 0);
	}
	return result;
}

float SHVector::findMinimum(int res = 32) const {
	float hExt = (float) M_PI / res, hInt = (2 * (float) M_PI)/(res*2);
	float minimum = std::numeric_limits<float>::infinity();

	for (int i=0; i<=res; ++i) {
		float theta = hExt*i;
		for (int j=0; j<=res*2; ++j) {
			float phi = hInt*j;
			minimum = std::min(minimum, eval(theta, phi));
		}
	}

	return minimum;
}

void SHVector::addOffset(float value) {
	operator()(0, 0) += 2 * value * (float) std::sqrt(M_PI);
}

float SHVector::eval(const Vector &v) const {
	float result = 0;
	float cosTheta = v.z, phi = std::atan2(v.y, v.x);
	if (phi < 0) phi += 2*M_PI;
	float *sinPhi = (float *) alloca(sizeof(float)*m_bands),
		  *cosPhi = (float *) alloca(sizeof(float)*m_bands);

	for (int m=0; m<m_bands; ++m) {
		sinPhi[m] = std::sin((m+1) * phi);
		cosPhi[m] = std::cos((m+1) * phi);
	}

	for (int l=0; l<m_bands; ++l) {
		for (int m=1; m<=l; ++m) {
			float L = legendreP(l, m, cosTheta) * normalization(l, m);
			result += operator()(l, -m) * SQRT_TWO * sinPhi[m-1] * L;
			result += operator()(l, m)  * SQRT_TWO * cosPhi[m-1] * L;
		}

		result += operator()(l, 0) * legendreP(l, 0, cosTheta) * normalization(l, 0);
	}
	return result;
}

float SHVector::evalAzimuthallyInvariant(float theta, float phi) const {
	float result = 0, cosTheta = std::cos(theta);
	for (int l=0; l<m_bands; ++l)
		result += operator()(l, 0) * legendreP(l, 0, cosTheta) * normalization(l, 0);
	return result;
}

float SHVector::evalAzimuthallyInvariant(const Vector &v) const {
	float result = 0, cosTheta = v.z;
	for (int l=0; l<m_bands; ++l)
		result += operator()(l, 0) * legendreP(l, 0, cosTheta) * normalization(l, 0);
	return result;
}

void SHVector::normalize() {
	float correction = 1/(2 * (float) std::sqrt(M_PI)*operator()(0,0));

	for (size_t i=0; i<(size_t) m_coeffs.size(); ++i)
		m_coeffs[i] *= correction;
}

void SHVector::convolve(const SHVector &kernel) {
	SAssert(kernel.getBands() == m_bands);

	for (int l=0; l<m_bands; ++l) {
		float alpha = std::sqrt(4 * (float) M_PI / (2*l + 1));
		for (int m=-l; m<=l; ++m)
			operator()(l, m) *= alpha * kernel(l, 0);
	}
}

Matrix3x3 SHVector::mu2() const {
	const float sqrt5o3 = std::sqrt((float) 5/ (float) 3);
	const float sqrto3 = std::sqrt((float) 1/ (float) 3);
	Matrix3x3 result;
	result.setZero();

	SAssert(m_bands > 0);
	result(0, 0) = result(1, 1) =
		result(2, 2) = sqrt5o3*operator()(0,0);

	if (m_bands >= 3) {
		result(0, 0) += -operator()(2,0)*sqrto3 + operator()(2,2);
		result(0, 1) = operator()(2,-2);
		result(0, 2) = -operator()(2, 1);
		result(1, 0) = operator()(2,-2);
		result(1, 1) += -operator()(2,0)*sqrto3 - operator()(2,2);
		result(1, 2) = -operator()(2,-1);
		result(2, 0) = -operator()(2, 1);
		result(2, 1) = -operator()(2,-1);
		result(2, 2) += 2*sqrto3*operator()(2,0);
	}

	return result * (2*std::sqrt((float) M_PI / 15));
}

std::string SHVector::toString() const {
	std::ostringstream oss;
	oss << "SHVector[bands=" << m_bands << ", {";
	int pos = 0;
	for (int i=0; i<m_bands; ++i) {
		oss << "{";
		for (int j=0, total=i*2+1; j<total; ++j) {
			oss << m_coeffs[pos++];
			if (j+1 < total)
				oss << ", ";
		}
		oss << "}";
		if (i+1 < m_bands)
			oss << ", ";
	}
	oss << "}]";
	return oss.str();
}

float SHVector::computeNormalization(int l, int m) {
	SAssert(m>=0);
	return std::sqrt(
			((2*l+1) * boost::math::factorial<float>(l-m))
		/    (4 * (float) M_PI * boost::math::factorial<float>(l+m)));
}

void SHVector::staticInitialization() {
	m_normalization = new float[SH_NORMTBL_SIZE*(SH_NORMTBL_SIZE+1)/2];
	for (int l=0; l<SH_NORMTBL_SIZE; ++l)
		for (int m=0; m<=l; ++m)
			m_normalization[l*(l+1)/2 + m] = computeNormalization(l, m);
}

void SHVector::staticShutdown() {
	delete[] m_normalization;
	m_normalization = NULL;
}

struct RotationBlockHelper {
	const SHRotation::Matrix &M1, &Mp;
	SHRotation::Matrix &Mn;
	int prevLevel, level;

	inline RotationBlockHelper(
		const SHRotation::Matrix &M1,
		const SHRotation::Matrix &Mp,
		SHRotation::Matrix &Mn)
		: M1(M1), Mp(Mp), Mn(Mn), prevLevel((int) Mp.rows()/2),
		level((int) Mp.rows()/2+1) { }

	inline float delta(int i, int j) const {
		return (i == j) ? (float) 1 : (float) 0;
	}

	inline float U(int l, int m, int n) const {
		return P(l, m, n, 0);
	}

	inline float V(int l, int m, int n) const {
		if (m == 0) {
			return P(l, 1, n, 1) + P(l, -1, n, -1);
		} else if (m > 0) {
			if (m == 1)
				return SQRT_TWO * P(l, 0, n, 1);
			else
				return P(l, m-1, n, 1) - P(l, -m+1, n, -1);
		} else {
			if (m == -1)
				return SQRT_TWO * P(l, 0, n, -1);
			else
				return P(l, -m-1, n, -1) + P(l, m+1, n, 1);
		}
	}

	inline float W(int l, int m, int n) const {
		if (m > 0) {
			return P(l, m+1, n, 1) + P(l, -m-1, n, -1);
		} else {
			return P(l, m-1, n, 1) - P(l, -m+1, n, -1);
		}
	}

	inline float u(int l, int m, int n) const {
		int denom = (std::abs(n) == l) ? (2*l*(2*l-1)) : ((l+n)*(l-n));
		return std::sqrt((float) ((l+m)*(l-m)) / (float) denom);
	}

	inline float v(int l, int m, int n) const {
		int denom = (std::abs(n) == l) ? (2*l*(2*l-1)) : ((l+n)*(l-n)), absM = std::abs(m);
		return .5f * (1-2*delta(m, 0)) * std::sqrt(
			(float) ((1+delta(m, 0)) * (l+absM-1)*(l+absM)) / (float) denom
		);
	}

	inline float w(int l, int m, int n) const {
		if (m == 0)
			return 0.0f;
		int absM = std::abs(m);
		int denom = (std::abs(n) < l) ? ((l+n)*(l-n)) : (2*l*(2*l-1));

		return -.5f * std::sqrt((float) ((l-absM-1)*(l-absM)) / (float) denom);
	}

	inline float P(int l, int m, int n, int i) const {
		if (std::abs(n) < l)
			return R(i, 0) * M(m, n);
		else if (n == l)
			return R(i, 1) * M(m, l-1) - R(i, -1) * M(m, -l+1);
		else if (n == -l)
			return R(i, 1) * M(m, -l+1) + R(i, -1) * M(m, l-1);
		else {
			SLog(EError, "Internal error!");
			return 0.0f;
		}
	}

	inline float R(int m, int n) const {
		return M1(m+1, n+1);
	}

	inline float M(int m, int n) const {
		return Mp(m+prevLevel, n+prevLevel);
	}

	void compute() {
		for (int m=-level; m<=level; ++m) {
			for (int n=-level; n<=level; ++n) {
				float uVal = u(level, m, n), vVal = v(level, m, n), wVal = w(level, m, n);
				Mn(m+level, n+level) =
					  (uVal != 0 ? (uVal * U(level, m, n)) : (float) 0)
					+ (vVal != 0 ? (vVal * V(level, m, n)) : (float) 0)
					+ (wVal != 0 ? (wVal * W(level, m, n)) : (float) 0);
			}
		}
	}
};

void SHVector::rotationBlock(
		const SHRotation::Matrix &M1,
		const SHRotation::Matrix &Mp,
		SHRotation::Matrix &Mn) {
	RotationBlockHelper rbh(M1, Mp, Mn);
	rbh.compute();
}

void SHVector::rotation(const Transform &t, SHRotation &rot) {
	rot.blocks[0](0, 0) = 1;
	if (rot.blocks.size() <= 1)
		return;

	const Matrix4x4 &trafo = t.getMatrix();
	rot.blocks[1](0, 0) =  trafo.m[1][1];
	rot.blocks[1](0, 1) = -trafo.m[2][1];
	rot.blocks[1](0, 2) =  trafo.m[0][1];
	rot.blocks[1](1, 0) = -trafo.m[1][2];
	rot.blocks[1](1, 1) =  trafo.m[2][2];
	rot.blocks[1](1, 2) = -trafo.m[0][2];
	rot.blocks[1](2, 0) =  trafo.m[1][0];
	rot.blocks[1](2, 1) = -trafo.m[2][0];
	rot.blocks[1](2, 2) =  trafo.m[0][0];

	if (rot.blocks.size() <= 2)
		return;

	for (size_t i=2; i<rot.blocks.size(); ++i)
		rotationBlock(rot.blocks[1], rot.blocks[i-1], rot.blocks[i]);
}

void SHRotation::operator()(const SHVector &source, SHVector &target) const {
	SAssert(source.getBands() == target.getBands());
	for (int l=0; l<source.getBands(); ++l) {
		const SHRotation::Matrix &M = blocks[l];
		for (int m1=-l; m1<=l; ++m1) {
			float result = 0;
			for (int m2=-l; m2<=l; ++m2)
				result += M(m1+l, m2+l)*source(l, m2);
			target(l, m1) = result;
		}
	}
}

SHSampler::SHSampler(int bands, int depth) : m_bands(bands), m_depth(depth) {
	m_phiMap = new float**[depth+1];
	m_legendreMap = new float**[depth+1];
	m_normalization = new float[m_bands*(m_bands+1)/2];
	m_dataSize = m_bands*(m_bands+1)/2;
	Assert(depth >= 1);

	for (int i=0; i<=depth; ++i) {
		int res = 1 << i;
		float zStep  = -2 / (float) res;
		float phiStep = 2 * (float) M_PI / (float) res;
		m_phiMap[i] = new float*[res];
		m_legendreMap[i] = new float*[res];

		for (int j=0; j<res; ++j) {
			m_phiMap[i][j] = phiIntegrals(phiStep*j, phiStep*(j+1));
			m_legendreMap[i][j] = legendreIntegrals(1+zStep*j, 1+zStep*(j+1));
		}
	}

	for (int l=0; l<m_bands; ++l) {
		for (int m=0; m<=l; ++m) {
			float normFactor = boost::math::tgamma_delta_ratio(
				(float) (l - m + 1), (float) (2 * m), boost::math::policies::policy<>());
			normFactor = std::sqrt(normFactor * (2 * l + 1) / (4 * (float) M_PI));
			if (m != 0)
				normFactor *= SQRT_TWO;
			m_normalization[I(l, m)] = normFactor;
		}
	}
}

std::string SHSampler::toString() const {
	std::ostringstream oss;
	oss << "SHSampler[bands=" << m_bands << ", depth=" << m_depth
		<< ", size=" << (m_dataSize*sizeof(double))/1024 << " KiB]";
	return oss.str();
}

float SHSampler::warp(const SHVector &f, Point2 &sample) const {
	int i = 0, j = 0;
	float integral = 0, integralRoot = integrate(0, 0, 0, f);

	for (int depth = 1; depth <= m_depth; ++depth) {
		/* Do not sample negative areas */
		float q00 = std::max(integrate(depth, i, j, f), (float) 0);
		float q10 = std::max(integrate(depth, i, j+1, f), (float) 0);
		float q01 = std::max(integrate(depth, i+1, j, f), (float) 0);
		float q11 = std::max(integrate(depth, i+1, j+1, f), (float) 0);

		float z1 = q00 + q10, z2 = q01 + q11, phi1, phi2;
		float zNorm = (float) 1 / (z1+z2);
		z1 *= zNorm; z2 *= zNorm;

		if (sample.x < z1) {
			sample.x /= z1;
			phi1 = q00; phi2 = q10;
			i <<= 1;
		} else {
			sample.x = (sample.x - z1) / z2;
			phi1 = q01; phi2 = q11;
			i = (i+1) << 1;
		}

		float phiNorm = (float) 1 / (phi1+phi2);
		float phi1Norm = phi1*phiNorm, phi2Norm = phi2*phiNorm;

		if (sample.y <= phi1Norm) {
			sample.y /= phi1Norm;
			j <<= 1;
			integral = phi1;
		} else {
			sample.y = (sample.y - phi1Norm) / phi2Norm;
			j = (j+1) << 1;
			integral = phi2;
		}
	}

	float zStep = -2 / (float) (1 << m_depth);
	float phiStep = 2 * (float) M_PI / (float) (1 << m_depth);
	i >>= 1; j >>= 1;

	float z = 1 + zStep * i + zStep * sample.x;
	sample.x = std::acos(z);
	sample.y = phiStep * j + phiStep * sample.y;

	/* PDF of sampling the mip-map bin */
	float pdfBin = integral/integralRoot;

	/* Density within the bin */
	float density = -1/(zStep*phiStep);

	return density*pdfBin;
}

SHSampler::~SHSampler() {
	for (int i=0; i<=m_depth; ++i) {
		int res = 1 << i;
		for (int j=0; j<res; ++j) {
			delete[] m_phiMap[i][j];
			delete[] m_legendreMap[i][j];
		}
		delete[] m_phiMap[i];
		delete[] m_legendreMap[i];
	}
	delete[] m_phiMap;
	delete[] m_legendreMap;
	delete[] m_normalization;
}

float SHSampler::integrate(int depth, int zBlock, int phiBlock, const SHVector &f) const {
	float result = 0;

	for (int l=0; l<m_bands; ++l) {
		for (int m=-l; m<=l; ++m) {
			float basisIntegral = m_normalization[I(l, std::abs(m))]
				* lookupIntegral(depth, zBlock, phiBlock, l, m);
			result += basisIntegral * f(l, m);
		}
	}
	return result;
}

float *SHSampler::phiIntegrals(float a, float b) {
	float *sinPhiA = new float[m_bands+1];
	float *sinPhiB = new float[m_bands+1];
	float *cosPhiA = new float[m_bands+1];
	float *cosPhiB = new float[m_bands+1];
	float *result = new float[2*m_bands+1];
	m_dataSize += 2*m_bands+1;

	cosPhiA[0] = 1; sinPhiA[0] = 0;
	cosPhiB[0] = 1; sinPhiB[0] = 0;
	cosPhiA[1] = std::cos(a);
	sinPhiA[1] = std::sin(a);
	cosPhiB[1] = std::cos(b);
	sinPhiB[1] = std::sin(b);

	for (int m=2; m<=m_bands; ++m) {
		sinPhiA[m] = 2*sinPhiA[m-1]*cosPhiA[1] - sinPhiA[m-2];
		sinPhiB[m] = 2*sinPhiB[m-1]*cosPhiB[1] - sinPhiB[m-2];

		cosPhiA[m] = 2*cosPhiA[m-1]*cosPhiA[1] - cosPhiA[m-2];
		cosPhiB[m] = 2*cosPhiB[m-1]*cosPhiB[1] - cosPhiB[m-2];
	}

	for (int m=-m_bands; m<=m_bands; ++m) {
		if (m == 0)
			result[P(m)] = b-a;
		else if (m > 0)
			result[P(m)] = (sinPhiB[m]-sinPhiA[m])/m;
		else
			result[P(m)] = (cosPhiB[-m]-cosPhiA[-m])/m;
	}

	delete[] sinPhiA;
	delete[] sinPhiB;
	delete[] cosPhiA;
	delete[] cosPhiB;
	return result;
}

float *SHSampler::legendreIntegrals(float a, float b) {
	float *P = new float[m_bands*(m_bands+1)/2];
	m_dataSize += m_bands*(m_bands+1)/2;

	P[I(0, 0)] = b-a;

	if (m_bands == 1)
		return P;

	float *Pa = new float[m_bands*(m_bands+1)/2];
	float *Pb = new float[m_bands*(m_bands+1)/2];

	for (int l=0; l<m_bands; ++l) {
		for (int m=0; m<=l; ++m) {
			Pa[I(l,m)] = legendreP(l, m, a);
			Pb[I(l,m)] = legendreP(l, m, b);
		}
	}

	P[I(1,0)] = (b*b - a*a)/2;
	P[I(1,1)] = .5f * (-b*std::sqrt(1-b*b) - std::asin(b) + a*std::sqrt(1-a*a) + std::asin(a));

	for (int l=2; l<m_bands; ++l) {
		for (int m=0; m<=l-2; ++m) {
			float ga = (2*l-1)*(1-a*a) * Pa[I(l-1,m)];
			float gb = (2*l-1)*(1-b*b) * Pb[I(l-1,m)];
			P[I(l, m)] = ((l-2)*(l-1+m)*P[I(l-2, m)]-gb+ga)/((l+1)*(l-m));
		}

		P[I(l, l-1)] = (2*l-1)/(float)(l+1) * ((1-a*a)*Pa[I(l-1, l-1)] - (1-b*b)*Pb[I(l-1, l-1)]);
		P[I(l, l)] = 1/(float)(l+1) * (l*(2*l-3)*(2*l-1) * P[I(l-2, l-2)] + b*Pb[I(l,l)] - a*Pa[I(l, l)]);
	}

	delete[] Pa;
	delete[] Pb;

	return P;
}

MTS_IMPLEMENT_CLASS(SHSampler, false, Object)
MTS_NAMESPACE_END
