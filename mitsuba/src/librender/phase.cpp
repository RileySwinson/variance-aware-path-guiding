#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>

MTS_NAMESPACE_BEGIN

std::string PhaseFunctionSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "PhaseFunctionSamplingRecord[" << endl
		<< "  mRec = " << indent(mRec.toString()) << "," << endl
		<< "  wi = " << wi.toString() << "," << endl
		<< "  wo = " << wo.toString() << "," << endl
		<< "  mode = " << mode << endl
		<< "]";
	return oss.str();
}

void PhaseFunction::configure() {
	m_type = 0;
}

float PhaseFunction::pdf(const PhaseFunctionSamplingRecord &pRec) const {
	return eval(pRec);
}

bool PhaseFunction::needsDirectionallyVaryingCoefficients() const {
	return false;
}

float PhaseFunction::sigmaDir(float cosTheta) const {
	Log(EError, "%s::sigmaDir(float) is not implemented (this is not "
		"an anisotropic medium!)", getClass()->getName().c_str());
	return 0.0f;
}

float PhaseFunction::sigmaDirMax() const {
	Log(EError, "%s::sigmaDirMax() is not implemented (this is not "
		"an anisotropic medium!)", getClass()->getName().c_str());
	return 0.0f;
}

float PhaseFunction::getMeanCosine() const {
	Log(EError, "%s::getMeanCosine() is not implemented!",
		getClass()->getName().c_str());
	return 0.0f;
}

MTS_IMPLEMENT_CLASS(PhaseFunction, true, ConfigurableObject)
MTS_NAMESPACE_END
