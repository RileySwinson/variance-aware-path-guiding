#include <mitsuba/ds/structures/gmm.h>

MTS_NAMESPACE_BEGIN

void GMM::construct(DSArguments& init_data)
{
    return;
}

void GMM::preprocess()
{
    return;
}

void GMM::store(std::vector<Sample>& samples)
{
    return;
}

void GMM::postprocess()
{
    return;
}

Sample GMM::sample(Point2& pos)
{
    Sample sample;
    sample.value = 0.0f;
    sample.pdf = 0.0f;
    sample.phi = 0.0f;
    sample.theta = 0.0f;

    return sample;
}

Float GMM::eval(Point2& pos)
{
    return 0;
}

void GMM::wipe()
{
    return;
}

DSType GMM::type()
{
    return DSType::DS_Invalid;
}

MTS_NAMESPACE_END