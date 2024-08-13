#include <mitsuba/ds/structures/unidir.h>

MTS_NAMESPACE_BEGIN

void Unidirectional::construct(DSArguments& init_data)
{
    this->sample_mode = init_data.mode;
}

void Unidirectional::preprocess()
{
    return;
}

void Unidirectional::store(std::vector<Sample>& samples)
{
    return; // We don't want to store anything for unguided unidirectional path tracing
}

void Unidirectional::postprocess()
{
    return;
}

Sample Unidirectional::sample(Point2& pos)
{
    // If the mode is set to cosine, we only want to return samples from the upper hemisphere.
    // Otherwise, return one over the the whole sphere.
    Vector dir = (this->sample_mode == Sample::Mode::Cosine)
        ? warp::squareToCosineHemisphere(pos)
        : warp::squareToUniformSphere(pos);

    Float theta = std::acos(dir.z);
    Float phi = std::atan2(dir.y, dir.x);
    if (phi < 0) phi += 2 * M_PI;

    //Float pdf = (this->sample_mode == Sample::Mode::Cosine)
    //    ? INV_PI * std::cos(theta)
    //    : INV_FOURPI;

    Sample sample;
    sample.value = 0;
    sample.pdf = 0;
    sample.phi = phi;
    sample.theta = theta;

    return sample;
}

Float Unidirectional::eval(Point2& pos)
{
    Point2 coords = Converter::uv_to_spherical(pos);

    //Float width = (0.5 * M_PI / 256) * (2 * M_PI / 1024);
    //Float width = (M_PI / 512) * (2 * M_PI / 1024);

    if (this->sample_mode == Sample::Mode::Cosine)
        return std::max((Float) 0, INV_PI * std::cos(coords.y) * std::sin(coords.y));

    return (INV_FOURPI * std::sin(coords.y));
}

void Unidirectional::wipe()
{
    return;
}

DSType Unidirectional::type()
{
    return DSType::DS_Unidirectional;
}

MTS_NAMESPACE_END