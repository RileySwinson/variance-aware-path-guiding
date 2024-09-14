#include <ds-compare/structures/vmfm.h>

MTS_NAMESPACE_BEGIN

void VMFM::construct(DSArguments& init_data)
{
    this->use_ruppert = init_data.use_ruppert;

    if (this->use_ruppert)
    {
        Properties props;
        props.setSize("vmmFactory.numInitialComponents", init_data.vmf_components);
        props.setFloat("vmmFactory.maxKappa", 32768.0f);
        props.setFloat("vmmFactory.rPriorWeight", 0.2f);
        props.setBoolean("parallaxCompensation", false);
        props.setBoolean("safetyMerge", false);

        this->factory_ruppert = VMMRuppertFactory(props);
    }
    else
    {
        VMMFactoryProperties vmm_native_props;
        vmm_native_props.numInitialComponents = init_data.vmf_components;
        vmm_native_props.maxKappa = 32768.0f;
        vmm_native_props.rPriorWeight = 0.2f;

        this->factory_native = VMMNativeFactory(vmm_native_props);
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
            VMMSample(Point3(), direction, sample.value, sample.pdf, INFINITY)
        );
    }
}

void VMFM::postprocess()
{
    if (!this->use_ruppert)
    {
        this->factory_native.fit(this->samples.begin(), this->samples.end(), this->vmm_native, true);
        return;
    }

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
    
    VMM4& vmm = (this->use_ruppert) 
        ? this->vmm_ruppert.distribution
        : this->vmm_native;

    return vmm.pdf(directional);
}

void VMFM::wipe()
{
    if (this->use_ruppert)
    {
        this->vmm_ruppert = VMMGuidingRegion();
    } 
    else
    {
        this->vmm_native = VMM4();
    }

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