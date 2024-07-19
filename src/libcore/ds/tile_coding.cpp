#include <mitsuba/ds/structures/tile_coding.h>

MTS_NAMESPACE_BEGIN

void TileCoding::construct(DSArguments& init_data)
{
    return;
}

void TileCoding::preprocess()
{
    return;
}

void TileCoding::store(Sample& sample)
{
    return;
}

void TileCoding::postprocess()
{
    return;
}

Sample TileCoding::sample(Point2& pos)
{
    Sample sample;
    sample.value = 0.0f;
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