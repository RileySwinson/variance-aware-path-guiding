#include <ds-compare/structures/vmfm.h>

MTS_NAMESPACE_BEGIN

void VMFM::construct(DSArguments& init_data)
{
    VMMFactoryProperties props; // TODO: Alter properties

    this->factory = VMMFactory(props);
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
        Vector3 direction( // Perhaps has to be inverted, tba
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
    this->factory.fit(this->samples.begin(), this->samples.end(), this->vmm, true);
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

    return this->vmm.pdf(directional);
}

void VMFM::wipe()
{
    this->vmm = VMM8();
    this->samples.clear();
}

DSType VMFM::type()
{
    return DSType::DS_VMFMixture;
}

MTS_NAMESPACE_END