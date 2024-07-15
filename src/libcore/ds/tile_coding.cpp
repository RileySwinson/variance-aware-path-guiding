#include <mitsuba/ds/structures/tile_coding.h>

MTS_NAMESPACE_BEGIN

DataStructure* TileCoding::construct(DSInitData& init_data)
{
    return this;
}

void TileCoding::store(Sample& sample)
{
    return;
}

Sample TileCoding::sample(Point2& pos)
{
    Sample sample;
    sample.luminance = 0.0f;
    sample.phi = 0.0f;
    sample.theta = 0.0f;

    return sample;
}

void TileCoding::wipe()
{
    return;
}

DSType TileCoding::type()
{
    return DSType::DS_TileCoding;
}

MTS_NAMESPACE_END