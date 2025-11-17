#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_VMFM_H_)
#define __DSCOMPARE_STRUCTURES_VMFM_H_

#include <pmm/DirectionalData.h>

#include <mitsuba/guiding/GuidingDistributionFactoryVMM.h>
#include <mitsuba/guiding/pathguidingmixturestatsfactory.h>
#include <mitsuba/guiding/GuidingFieldFactory.h>
#include <mitsuba/guiding/pathguidingmixturestats.h>

#include <ds-compare/ds.h>
#include <boost/variant.hpp>

MTS_NAMESPACE_BEGIN

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

struct NativeStrategy {
    VMMNativeFactory factory;
    VMM4 vmm;
};

struct RuppertStrategy {
    VMMRuppertFactory factory;
    VMMGuidingRegion vmm;
};

struct VMMStrategy {
    using VMMVariant = boost::variant<NativeStrategy, RuppertStrategy>;
    VMMStrategy() : m_strategy(NativeStrategy()) { };
    VMMStrategy(VMMVariant strategy) : m_strategy(strategy) { };

    VMM4& vmm() // Do not trust the linter. Do not make this const.
    {
        struct VMM4Visitor : public boost::static_visitor<VMM4&>
        {
            VMM4& operator()(NativeStrategy& strategy) const { return strategy.vmm; }
            VMM4& operator()(RuppertStrategy& strategy) const { return strategy.vmm.distribution; }
        };

        return boost::apply_visitor(VMM4Visitor(), this->m_strategy);
    }

    VMMVariant& get() { return this->m_strategy; }
    const VMMVariant& get() const { return this->m_strategy; }
private:
    VMMVariant m_strategy;
};

struct FitVisitor : public boost::static_visitor<>
{
    std::vector<VMMSample>& samples_ref;
    mutable bool first_fit = true;
    FitVisitor(std::vector<VMMSample>& samples) : samples_ref(samples) {}

    void operator()(NativeStrategy& strategy) const {
        if (this->first_fit)
        {
            strategy.factory.fit(this->samples_ref.begin(), this->samples_ref.end(), strategy.vmm, true);
            this->first_fit = false;
            return;
        }

        strategy.factory.updateFit(this->samples_ref.begin(), this->samples_ref.end(), strategy.vmm);
    }
    void operator()(RuppertStrategy& strategy) const {
        if (this->first_fit)
        {
            strategy.factory.fit(strategy.vmm, this->samples_ref);
            this->first_fit = false;
            return;
        }
        
        strategy.factory.updateFit(strategy.vmm, this->samples_ref);
    }
};

struct WipeVisitor : public boost::static_visitor<>
{
    void operator()(NativeStrategy& strategy) const {
        strategy.vmm = VMM4();
    }
    void operator()(RuppertStrategy& strategy) const {
        strategy.vmm = VMMGuidingRegion();
    }
};

struct MTS_EXPORT_CORE VMFM : public DataStructure {
    ~VMFM() { }

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
    VMMStrategy strategy;
    std::vector<VMMSample> curr_samples;
    float sample_memory = 0.0f; // Only used for memory estimation.
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_VMFM_H_ */