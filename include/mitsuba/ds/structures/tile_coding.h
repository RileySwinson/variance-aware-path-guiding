#pragma once
#if !defined(__MITSUBA_DS_STRUCTURES_TILE_CODING_H_)
#define __MITSUBA_DS_STRUCTURES_TILE_CODING_H_

#include <mitsuba/ds/ds.h>
#include <Eigen/Core>

MTS_NAMESPACE_BEGIN

struct MTS_EXPORT_CORE TileCoding : public DataStructure {
    ~TileCoding() { }

    void construct(DSArguments& init_data) override;

    void preprocess() override;

    void store(Sample& sample) override;

    void postprocess() override;

    Sample sample(Point2& pos) override;
    
    void wipe() override;

    DSType type() override;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_DS_STRUCTURES_TILE_CODING_H_ */