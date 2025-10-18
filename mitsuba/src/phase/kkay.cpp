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

#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{kkay}{Kajiya-Kay phase function}
 * This plugin implements the Kajiya-Kay \cite{Kajiya1989Rendering}
 * phase function for volumetric rendering of fibers, e.g.
 * hair or cloth.
 *
 * The function is normalized so that it has no energy loss when
 * \code{ks}=1 and illumination arrives perpendicularly to the surface.
 */
class KajiyaKayPhaseFunction : public PhaseFunction {
public:
	KajiyaKayPhaseFunction(const Properties &props)
		: PhaseFunction(props) {
		m_ks = props.getfloat("ks", .4f);
		m_kd = props.getfloat("kd", .2f);
		m_exponent = props.getfloat("exponent", 4.0f);
		if (m_kd + m_ks > 1.0f)
			Log(EWarn, "Energy conservation is violated!");
	}

	KajiyaKayPhaseFunction(Stream *stream, InstanceManager *manager)
		: PhaseFunction(stream, manager) {
		m_ks = stream->readfloat();
		m_kd = stream->readfloat();
		m_exponent = stream->readfloat();
		configure();
	}

	virtual ~KajiyaKayPhaseFunction() { }

	void configure() {
		/* Compute the normalization for perpendicular illumination
		   using Simpson quadrature */
		int nParts = 1000;
		float stepSize = M_PI / nParts, m=4, theta = stepSize;

		m_normalization = 0; /* 0 at the endpoints */
		for (int i=1; i<nParts; ++i) {
			float value = std::pow(std::cos(theta - M_PI/2), m_exponent)
				* std::sin(theta);
			m_normalization += value * m;
			theta += stepSize;
			m = 6-m;
		}

		m_normalization = 1/(m_normalization * stepSize/3 * 2 * M_PI);
		m_type = EAnisotropic;
		Log(EDebug, "Kajiya-kay normalization factor = %f", m_normalization);
	}


	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
		stream->writefloat(m_ks);
		stream->writefloat(m_kd);
		stream->writefloat(m_exponent);
	}

	float sample(PhaseFunctionSamplingRecord &pRec,
			Sampler *sampler) const {
		pRec.wo = warp::squareToUniformSphere(sampler->next2D());
		return eval(pRec) * (4 * M_PI);
	}

	float sample(PhaseFunctionSamplingRecord &pRec,
			float &pdf, Sampler *sampler) const {
		pRec.wo = warp::squareToUniformSphere(sampler->next2D());
		pdf = warp::squareToUniformSpherePdf();
		return eval(pRec) * (4 * M_PI);
	}

	float pdf(const PhaseFunctionSamplingRecord &pRec) const {
		return warp::squareToUniformSpherePdf();
	}

	float eval(const PhaseFunctionSamplingRecord &pRec) const {
		if (pRec.mRec.orientation.length() == 0)
			return m_kd / (4*M_PI);

		Frame frame(normalize(pRec.mRec.orientation));
		Vector reflectedLocal = frame.toLocal(pRec.wo);

		reflectedLocal.z = -dot(pRec.wi, frame.n);
		float a = std::sqrt((1-reflectedLocal.z*reflectedLocal.z) /
			(reflectedLocal.x*reflectedLocal.x + reflectedLocal.y*reflectedLocal.y));
		reflectedLocal.y *= a;
		reflectedLocal.x *= a;
		Vector R = frame.toWorld(reflectedLocal);

		return std::pow(std::max((float) 0, dot(R, pRec.wo)), m_exponent)
			* m_normalization * m_ks + m_kd / (4*M_PI);
	}

	std::string toString() const {
		return "KajiyaKayPhaseFunction[]";
	}

	MTS_DECLARE_CLASS()
private:
	float m_ks, m_kd, m_exponent, m_normalization;
};


MTS_IMPLEMENT_CLASS_S(KajiyaKayPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(KajiyaKayPhaseFunction, "Kajiya-Kay phase function");
MTS_NAMESPACE_END
