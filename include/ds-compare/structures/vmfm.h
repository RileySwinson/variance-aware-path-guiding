#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_VMFM_H_)
#define __DSCOMPARE_STRUCTURES_VMFM_H_

#include <pmm/DirectionalData.h>

#include <mitsuba/guiding/GuidingDistributionFactoryVMM.h>
#include <mitsuba/guiding/pathguidingmixturestatsfactory.h>
#include <mitsuba/guiding/GuidingFieldFactory.h>
#include <mitsuba/guiding/pathguidingmixturestats.h>

#include <ds-compare/ds.h>

MTS_NAMESPACE_BEGIN

enum VMMMode {
    Native,
    Ruppert
};

// I mean, I get it. Ruppert et al. wanted this code to be as flexible as possible,
// which is why everyone and their mom got a template. However, this makes dynamic
// usage *very* painful. Here we go...
using SSEScalar = ::lightpmm::Scalar4;
using VMFKernel = ::lightpmm::VMFKernel<SSEScalar>;
using VMM4 = ::lightpmm::ParametricMixtureModel<VMFKernel, 4>;
using VMMStatistics = Guiding::PathGuidingMixtureStats<VMM4>;

using VMMGuidingRegion = Guiding::GuidingRegion<VMM4, VMMStatistics>;
using VMMNativeFactory = ::lightpmm::VMMFactory<VMM4>;
using VMMRuppertFactory = Guiding::GuidingFieldFactory<VMM4, VMMStatistics>;

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
    VMMMode mode;

    // Native factory + vmm
    VMMNativeFactory factory_native;
    VMM4 vmm_native;

    // Ruppert factory + vmm
    VMMRuppertFactory factory_ruppert;
    VMMGuidingRegion vmm_ruppert;

    std::vector<VMMSample> samples;
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_VMFM_H_ */