#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_VMFM_H_)
#define __DSCOMPARE_STRUCTURES_VMFM_H_

#include <pmm/DirectionalData.h>
#include <pmm/VMMFactory.h>

#include <ds-compare/ds.h>

MTS_NAMESPACE_BEGIN

// I mean, I get it. Ruppert et al. wanted this code to be as flexible as possible,
// which is why everyone and their mom got a template. However, this makes dynamic
// usage *very* painful. Here we go...
using Scalar8 = ::lightpmm::Scalar8;
using VMFKernel = ::lightpmm::VMFKernel<Scalar8>;
using VMM8 = ::lightpmm::ParametricMixtureModel<VMFKernel, 8>;
using VMMFactory = ::lightpmm::VMMFactory<VMM8>;

using VMMFactoryProperties = ::lightpmm::VMMFactoryProperties;
using VMMSample = ::lightpmm::DirectionalData;

struct MTS_EXPORT_CORE VMFM : public DataStructure {
    ~VMFM() { }

    void construct(DSArguments& init_data) override;

    void preprocess() override;

    void store(std::vector<Sample>& samples) override;

    void postprocess() override;

    Sample sample(Point2& pos) override;

    Float eval(Point2& pos) override;

    void wipe() override;

    DSType type() override;

    std::string name() override;

    int memory() override;

private:
    VMMFactory factory;
    VMM8 vmm;
    std::vector<VMMSample> samples;
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_VMFM_H_ */