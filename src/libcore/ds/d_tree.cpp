#include <mitsuba/ds/structures/d_tree.h>

MTS_NAMESPACE_BEGIN

DataStructure* DirectionalTree::construct(DSInitData& init_data)
{
    return this;
}

void DirectionalTree::store(Sample& sample)
{
    return;
}

Sample DirectionalTree::sample(Point2& pos)
{
    Sample sample;
    sample.luminance = 0.0f;
    sample.phi = 0.0f;
    sample.theta = 0.0f;

    return sample;
}

void DirectionalTree::wipe()
{
    return;
}

DSType DirectionalTree::type()
{
    return DSType::DS_Invalid;
}

MTS_NAMESPACE_END