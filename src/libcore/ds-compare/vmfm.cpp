#include <ds-compare/structures/vmfm.h>

MTS_NAMESPACE_BEGIN

void VMFM::construct(DSArguments& init_data)
{
    if (init_data.vmf.use_ruppert)
    {
        Properties props;
        props.setSize("vmmFactory.numInitialComponents", init_data.vmf.components);
        props.setFloat("vmmFactory.maxKappa", 32768.0f);
        props.setFloat("vmmFactory.rPriorWeight", 0.2f);
        props.setBoolean("parallaxCompensation", false);
        props.setBoolean("safetyMerge", false);

        this->strategy = VMMStrategy(RuppertStrategy{ VMMRuppertFactory(props), VMMGuidingRegion() });
    }
    else
    {
        VMMFactoryProperties vmm_native_props;
        vmm_native_props.numInitialComponents = init_data.vmf.components;
        vmm_native_props.maxKappa = 32768.0f;
        vmm_native_props.rPriorWeight = 0.2f;

        this->strategy = VMMStrategy(NativeStrategy{ VMMNativeFactory(vmm_native_props), VMM4() });
    }
}

void VMFM::preprocess()
{
    return;
}

void VMFM::store(std::vector<Sample>& input_samples)
{
    std::vector<VMMSample> samples;

    for (const auto& sample : input_samples)
    {
        Vector3 direction(
            std::sin(sample.theta) * std::cos(sample.phi),
            std::sin(sample.theta) * std::sin(sample.phi),
            std::cos(sample.theta)
        );

        samples.emplace_back(
            VMMSample(Point3(), direction, sample.value, sample.pdf, Epsilon)
        );
    }

    boost::apply_visitor(FitVisitor(samples), this->strategy.get());
}

void VMFM::postprocess()
{
    return;
}

Sample VMFM::sample(Point2& pos)
{
    VMM4& vmm = this->strategy.vmm();

    Vector3 dir = vmm.sample(pos);
    float pdf = vmm.pdf(dir);

    Float theta = std::acos(dir.z);
    Float phi = std::atan2(dir.y, dir.x);
    if (phi < 0) phi += 2 * M_PI;

    Sample sample;
    sample.value = 0.0f;
    sample.pdf = pdf;
    sample.phi = boost::algorithm::clamp(phi, 0, 2 * M_PI - Epsilon);
    sample.theta = boost::algorithm::clamp(theta, 0, M_PI - Epsilon);

    return sample;
}

Float VMFM::eval(Point2& pos)
{
    Point2 spherical = Converter::uv_to_spherical(pos);
    Vector3 directional(
        std::sin(spherical.y) * std::cos(spherical.x),
        std::sin(spherical.y) * std::sin(spherical.x),
        std::cos(spherical.y)
    );
    
    VMM4& vmm = this->strategy.vmm();
    return vmm.pdf(directional);
}

void VMFM::wipe()
{
    boost::apply_visitor(WipeVisitor(), this->strategy.get());
}

std::string VMFM::name()
{
    return "Von-Mises-Fisher Mixture";
}

DSType VMFM::type()
{
    return DSType::DS_VMFMixture;
}

int VMFM::memory()
{
    // It completely suffices to call sizeof, as there is no dynamic storage or
    // similar used in this vMFM implementation. This obviously only serves as a
    // lower bound, as there still may be platform-dependent dynamic allocations
    // happening, but this generally applies to all memory approximations.

    return sizeof(this->strategy.vmm());
}

MTS_NAMESPACE_END