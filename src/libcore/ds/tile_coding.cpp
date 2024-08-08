#include <mitsuba/ds/structures/tile_coding.h>

MTS_NAMESPACE_BEGIN

void TileCoding::construct(DSArguments& init_data)
{
    for (auto tiling : this->m_tilings)
    {
        tiling.resize(this->m_tiles * this->m_tiles);
    }
}

void TileCoding::preprocess()
{
    return;
}

void TileCoding::store(std::vector<Sample>& samples)
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

Float TileCoding::eval(Point2& pos)
{
    return 0;
}

void TileCoding::wipe()
{
    return;
}

DSType TileCoding::type()
{
    return DSType::DS_Invalid;
}

MTS_NAMESPACE_END