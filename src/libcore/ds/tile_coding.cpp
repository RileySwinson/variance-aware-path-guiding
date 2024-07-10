#include <mitsuba/ds/structures/tile_coding.h>

MTS_NAMESPACE_BEGIN

DataStructure* TileCoding::construct()
{
    return this;
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