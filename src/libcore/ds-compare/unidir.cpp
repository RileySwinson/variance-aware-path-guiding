#include <ds-compare/structures/unidir.h>

MTS_NAMESPACE_BEGIN

void Unidirectional::construct_impl(DSArguments& init_data)
{
    this->sample_mode = init_data.comparer.mode;
}

void Unidirectional::preprocess_impl()
{
    return;
}

void Unidirectional::store_impl(Sample& sample)
{
    return; // We don't want to store anything for unguided unidirectional path tracing
}

void Unidirectional::postprocess_impl(bool last_iteration)
{
    return;
}

Sample Unidirectional::sample_impl(Point2& pos)
{
    // If the mode is set to cosine, we only want to return samples from the upper hemisphere.
    // Otherwise, return one over the the whole sphere.
    Vector dir = (this->sample_mode == Sample::Mode::Cosine)
        ? warp::squareToCosineHemisphere(pos)
        : warp::squareToUniformSphere(pos);

    Float theta = std::acos(dir.z);
    Float phi = std::atan2(dir.y, dir.x);
    while (phi < 0) phi += 2 * M_PI;

    Point2 spherical(phi, theta);

    Sample sample;
    sample.value = 0;
    sample.pdf = pdf(spherical);
    sample.phi = phi;
    sample.theta = theta;

    return sample;
}

Float Unidirectional::eval_impl(Point2& pos)
{
    Point2 coords = Converter::uv_to_spherical(pos);
    return pdf(coords);
}

void Unidirectional::wipe_impl()
{
    return;
}

DSType Unidirectional::type_impl()
{
    return DSType::DS_Unidirectional;
}

std::string Unidirectional::name_impl()
{
    return "Unidirectional";
}

Float Unidirectional::pdf(Point2& coords) const
{
    if (this->sample_mode == Sample::Mode::Cosine)
        return std::max((Float) 0, INV_PI * std::cos(coords.y));

    return INV_FOURPI;
}

int Unidirectional::memory_impl()
{
    return 0; // In a real path tracing setting no data structure is used, hence we return 0.
}

MTS_NAMESPACE_END