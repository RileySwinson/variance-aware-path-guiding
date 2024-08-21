#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_NONE_H_)
#define __DSCOMPARE_STRUCTURES_NONE_H_

#include <ds-compare/ds.h>
#include <Eigen/Core>

MTS_NAMESPACE_BEGIN

struct MTS_EXPORT_CORE Unidirectional : public DataStructure {
    ~Unidirectional() { }

    void construct(DSArguments& init_data) override;

    void preprocess() override;

    void store(std::vector<Sample>& samples) override;

    void postprocess() override;

    Sample sample(Point2& pos) override;

    Float eval(Point2& pos) override;

    void wipe() override;

    DSType type() override;

    std::string name() override;

private:
    Sample::Mode sample_mode;
    ref<Random> random = new Random();

    /// Takes a position in spherical coordinates [phi, theta] and returns the pdf at that position.
    Float pdf(Point2& pos);
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_NONE_H_ */