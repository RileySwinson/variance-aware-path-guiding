#pragma once
#if !defined(__MITSUBA_DS_STRUCTURES_TILE_CODING_H_)
#define __MITSUBA_DS_STRUCTURES_TILE_CODING_H_

#include <mitsuba/ds/ds.h>
#include <Eigen/Core>

MTS_NAMESPACE_BEGIN

struct MTS_EXPORT_CORE TileCoding : public DataStructure {
    ~TileCoding() { }

    DataStructure* construct() override;

    void wipe() override;

    DSType type() override;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_DS_STRUCTURES_TILE_CODING_H_ */