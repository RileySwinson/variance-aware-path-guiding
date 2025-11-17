#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_NONE_H_)
#define __DSCOMPARE_STRUCTURES_NONE_H_

#include <ds-compare/ds.h>

MTS_NAMESPACE_BEGIN

struct MTS_EXPORT_CORE Unidirectional : public DataStructure {
    ~Unidirectional() { }

    void construct_impl(DSArguments& init_data) override;

    void preprocess_impl() override;

    void store_impl(Sample& sample) override;

    void postprocess_impl(bool last_iteration = false) override;

    Sample sample_impl(Point2& pos) override;

    Float eval_impl(Point2& pos) override;

    void wipe_impl() override;

    DSType type_impl() override;

    std::string name_impl() override;

    int memory_impl() override;

private:
    Sample::Mode sample_mode;
    ref<Random> random = new Random();

    /// Takes a position in spherical coordinates [phi, theta] and returns the pdf at that position.
    Float pdf(Point2& pos) const;
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_NONE_H_ */