#include <ds-compare/structures/tile_coding_dyn.h>

MTS_NAMESPACE_BEGIN

void DynamicTileCoding::construct(DSArguments& init_data)
{
    // TODO
    this->m_mode = init_data.mode;
}

void DynamicTileCoding::preprocess()
{
    // TODO
    return;
}

void DynamicTileCoding::store(std::vector<Sample>& samples)
{
    // TODO
    return;
}

void DynamicTileCoding::postprocess()
{
    // TODO
    return;
}

Sample DynamicTileCoding::sample(Point2& pos)
{
    // TODO
    Sample sample;
    sample.value = 0;
    sample.pdf = 0;
    sample.phi = 0;
    sample.theta = 0;

    return sample;
}

Float DynamicTileCoding::eval(Point2& pos)
{
    return pdf(pos);
}

void DynamicTileCoding::wipe()
{
    return;
}

DSType DynamicTileCoding::type()
{
    return DSType::DS_DynamicTileCoding;
}

std::string DynamicTileCoding::name()
{
    return "Dynamic Tile Coding";
}

Float DynamicTileCoding::pdf(Point2& coords)
{
    // TODO
    return 0;
}

int DynamicTileCoding::memory()
{
    // TODO
    return 0;
}

MTS_NAMESPACE_END