#pragma once
#if !defined(__MITSUBA_DS_STRUCTURES_D_TREE_H_)
#define __MITSUBA_DS_STRUCTURES_D_TREE_H_

#include <mitsuba/ds/ds.h>
#include <Eigen/Core>

MTS_NAMESPACE_BEGIN

struct MTS_EXPORT_CORE DirectionalTree : public DataStructure {
    ~DirectionalTree() { }

    void construct(DSInitData& init_data) override;

    void store(Sample& sample) override;

    Sample sample(Point2& pos) override;
    
    void wipe() override;

    DSType type() override;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_DS_STRUCTURES_D_TREE_H_ */