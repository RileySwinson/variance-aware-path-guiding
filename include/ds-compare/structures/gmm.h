#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_GMM_H_)
#define __DSCOMPARE_STRUCTURES_GMM_H_

#include <ds-compare/ds.h>
#include <Eigen/Core>

MTS_NAMESPACE_BEGIN

struct MTS_EXPORT_CORE GMM : public DataStructure {
    ~GMM() { }

    void construct(DSArguments& init_data) override;

    void preprocess() override;

    void store(std::vector<Sample>& samples) override;

    void postprocess() override;

    Sample sample(Point2& pos) override;

    Float eval(Point2& pos) override;

    void wipe() override;

    DSType type() override;
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_GMM_H_ */