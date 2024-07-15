#include <mitsuba/ds/structures/gmm.h>

MTS_NAMESPACE_BEGIN

DataStructure* GMM::construct(DSInitData& init_data)
{
    return this;
}

void GMM::store(Sample& sample)
{
    return;
}

Sample GMM::sample(Point2& pos)
{
    Sample sample;
    sample.luminance = 0.0f;
    sample.phi = 0.0f;
    sample.theta = 0.0f;

    return sample;
}

void GMM::wipe()
{
    return;
}

DSType GMM::type()
{
    return DSType::DS_GaussianMixture;
}

MTS_NAMESPACE_END