#include <ds-compare/structures/vmfm.h>

MTS_NAMESPACE_BEGIN

void VMFM::construct(DSArguments& init_data)
{
    this->mode = VMMMode::Native;

    if (this->mode == VMMMode::Native)
    {
        VMMFactoryProperties vmm_native_props;
        vmm_native_props.numInitialComponents = init_data.vmf_components;
        vmm_native_props.maxKappa = 32768.0f;
        vmm_native_props.rPriorWeight = 0.2f;

        this->factory_native = VMMNativeFactory(vmm_native_props);
    }

    if (this->mode == VMMMode::Ruppert)
    {
        Properties vmm_ruppert_props;
        vmm_ruppert_props.setBoolean("parallaxCompensation", false);

        this->factory_ruppert = VMMRuppertFactory(vmm_ruppert_props);
    }
    
    samples.reserve(init_data.samples_learning);
}

void VMFM::preprocess()
{
    return;
}

void VMFM::store(std::vector<Sample>& input_samples)
{
    for (const auto& sample : input_samples)
    {
        Vector3 direction(
            std::sin(sample.theta) * std::cos(sample.phi),
            std::sin(sample.theta) * std::sin(sample.phi),
            std::cos(sample.theta)
        );

        this->samples.emplace_back(
            VMMSample(Point3(), direction, sample.value, sample.pdf, 0.0f)
        );
    }
}

void VMFM::postprocess()
{
    if (this->mode == VMMMode::Native)
        this->factory_native.fit(this->samples.begin(), this->samples.end(), this->vmm_native, true);

    if (this->mode == VMMMode::Ruppert)
        this->factory_ruppert.fit(this->vmm_ruppert, this->samples);
}

Sample VMFM::sample(Point2& pos)
{
    Sample sample;
    sample.value = 0.0f;
    sample.pdf = 0.0f;
    sample.phi = 0.0f;
    sample.theta = 0.0f;

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
    
    VMM4& vmm = (this->mode == VMMMode::Native) 
        ? this->vmm_native 
        : this->vmm_ruppert.distribution;

    return vmm.pdf(directional);
}

void VMFM::wipe()
{
    if (this->mode == VMMMode::Native)
        this->vmm_native = VMM4();

    if (this->mode == VMMMode::Ruppert)
        this->vmm_ruppert = VMMGuidingRegion();

    this->samples.clear();
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
    // TODO
    return 0;
}

MTS_NAMESPACE_END