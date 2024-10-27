#pragma once

#if !defined(__DSCOMPARE_UTIL_RANDOM_H_)
#define __DSCOMPARE_UTIL_RANDOM_H_

#include <ds-compare/util/_definitions.h>

MTS_NAMESPACE_BEGIN

struct DS_COMPARE RandomGen {
    Float next1D() {
        return m_random->nextFloat();
    }

    Point2 next2D() {
        Float v1 = m_random->nextFloat();
        Float v2 = m_random->nextFloat();
        return Point2(v1, v2);
    }

private:
    ref<Random> m_random = new Random();
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_UTIL_RANDOM_H_ */