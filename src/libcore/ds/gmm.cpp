#include <mitsuba/ds/structures/gmm.h>

MTS_NAMESPACE_BEGIN

DataStructure* GMM::construct()
{
    return this;
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