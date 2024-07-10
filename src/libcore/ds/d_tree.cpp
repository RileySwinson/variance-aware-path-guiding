#include <mitsuba/ds/structures/d_tree.h>

MTS_NAMESPACE_BEGIN

DataStructure* DirectionalTree::construct()
{
    return this;
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