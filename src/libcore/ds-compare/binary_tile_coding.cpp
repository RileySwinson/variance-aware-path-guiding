#include <ds-compare/structures/binary_tile_coding.h>
#include <iomanip>

MTS_NAMESPACE_BEGIN

/* ========== */
/* BinaryTile */
/* ========== */

bool BinaryTile::is_leaf() const
{
    return (this->idx_first == -1) && (this->idx_second == -1);
}

/// Checks if the current tile fulfills all criteria for splitting.
/// Make sure that this tile isn't a leaf by checking against !is_leaf() before!
bool BinaryTile::is_splittable() const
{
    // Are there enough samples in this cell?
    if (this->meta_data.sample_count < BinaryTileCoding::MIN_SAMPLES)
    {
        return false;
    }

    return (std::abs(covar(HORIZONTAL)) > BinaryTileCoding::SUBDIV_THRESHOLD 
        || std::abs(covar(VERTICAL)) > BinaryTileCoding::SUBDIV_THRESHOLD);
}

SplitDirection BinaryTile::split_direction() const
{
    return static_cast<SplitDirection>(this->meta_data.split);
}

void BinaryTile::update_covariance(Sample& sample)
{
    auto samples = this->meta_data.sample_count;
    float value_mean = (samples > 0)
        ? (this->value / samples)
        : 0;

    auto n = samples + 1;

    float dx = sample.value - value_mean;
    auto dy_x = [&sample, this]() { return sample.phi - this->sample_mean.x; };
    auto dy_y = [&sample, this]() { return sample.theta - this->sample_mean.y; };

    // Update x covariance
    this->sample_mean.x += dy_x() / n;
    this->cov.x += dx * dy_x();

    // Update y covariance
    this->sample_mean.y += dy_y() / n;
    this->cov.y += dx * dy_y();
}

void BinaryTile::update_value(Sample& sample)
{
    this->value += sample.value;
    this->meta_data.sample_count++;
}

float BinaryTile::covar(SplitDirection dir) const
{
    float c = (dir == HORIZONTAL) ? this->cov.x : this->cov.y;
    return c / (this->meta_data.sample_count - 1);
}

/* ============ */
/* BinaryTiling */
/* ============ */

void BinaryTiling::insert(Sample& sample, Point2i& tile_dims)
{
    // Find initial tile
    Point2 uv = Converter::spherical_to_uv(Point2(sample.phi, sample.theta));

    Float x_pos = this->factors.x * (uv.x - this->start_vals.x);
    Float y_pos = this->factors.y * (uv.y - this->start_vals.y);

    Point2i index(
        x_pos * tile_dims.x,
        y_pos * tile_dims.y
    );

    int i = (index.y * tile_dims.x) + index.x;
    BinaryTile& curr_tile = this->tiles.at(i);

    // Calculate boundary of current tile
    Point2 x_bounds(index.x / (Float) tile_dims.x, (index.x + 1) / (Float) tile_dims.x);
    Point2 y_bounds(index.y / (Float) tile_dims.y, (index.y + 1) / (Float) tile_dims.y);

    // Iterate through tree if necessary
    while (!curr_tile.is_leaf())
    {
        SplitDirection split_dir = curr_tile.split_direction();

        Point2& bounds = (split_dir == VERTICAL) ? x_bounds : y_bounds;
        Float& pos = (split_dir == VERTICAL) ? x_pos : y_pos;
        Float split = (bounds.x + bounds.y) * 0.5;

        if (pos < split) // we're in the left or upper subtree
        {
            bounds.y = split;
            curr_tile = this->tiles.at(curr_tile.idx_first);
        }
        else // we're in the right or lower subtree
        {
            bounds.x = split;
            curr_tile = this->tiles.at(curr_tile.idx_second);
        }
    }

    Sample updater;
    updater.value = sample.value;
    updater.phi = x_pos;
    updater.theta = y_pos;

    // Store value & update covariance
    curr_tile.update_covariance(updater);
    curr_tile.update_value(updater);
    
    // Split if necessary
    if (!curr_tile.is_splittable()) return;

    this->tiles.push_back(BinaryTile());
    this->tiles.push_back(BinaryTile());

    curr_tile.idx_first = this->tiles.size();
    curr_tile.idx_second = curr_tile.idx_first + 1;
    
    curr_tile.meta_data.split = (std::abs(curr_tile.covar(VERTICAL)) > BinaryTileCoding::SUBDIV_THRESHOLD) 
        ? VERTICAL 
        : HORIZONTAL;
}

/* ================ */
/* BinaryTileCoding */
/* ================ */

void BinaryTileCoding::construct(DSArguments& init_data)
{
    this->tile_dims = Point2i(init_data.tiles_x, init_data.tiles_y);

    this->tilings.resize(init_data.tilings);
    for (auto& tiling : this->tilings)
    {
        tiling.tiles.resize(this->tile_dims.x * this->tile_dims.y);
    }
}

void BinaryTileCoding::preprocess()
{
    // Calculate tiling offset
    Float tile_width = 1.0 / this->tile_dims.x;
    Float tile_height = 1.0 / this->tile_dims.y;
    auto num_tilings = this->tilings.size();

    Float overhead = 0.0;
    if (num_tilings > 1)
    {
        overhead = 1.0 / num_tilings;
    }

    Point2 offset(tile_width * overhead, tile_height * overhead);

    // Store start values and mapping factors in each tiling
    for (int i = 0; i < num_tilings; ++i)
    {
        BinaryTiling& tiling = this->tilings.at(i);

        Point2 x_range(0 - (i * offset.x), 1 + ((num_tilings - 1 - i) * offset.x));
        Point2 y_range(0 - ((num_tilings - 1 - i) * offset.y), 1 + (i * offset.y));

        Float x_factor = 1.0 / (x_range.y - x_range.x);
        Float y_factor = 1.0 / (y_range.y - y_range.x);

        tiling.start_vals = Point2(x_range.x, y_range.x);
        tiling.factors = Point2(x_factor, y_factor);
    }
}

void BinaryTileCoding::store(std::vector<Sample>& samples)
{
    for (auto& sample : samples)
    {
        for (auto& tiling : this->tilings)
        {
            tiling.insert(sample, this->tile_dims);
        }
    }
}

void BinaryTileCoding::postprocess()
{
    // TODO: Create stuff for eval & sampling
    return;
}

Sample BinaryTileCoding::sample(Point2& pos)
{
    Sample sample;
    return sample;
}

Float BinaryTileCoding::eval(Point2& pos)
{
    return 0.0f;
}

void BinaryTileCoding::wipe()
{
    return;
}

DSType BinaryTileCoding::type()
{
    return DSType::DS_BinaryTileCoding;
}

std::string BinaryTileCoding::name()
{
    return "Binary Tile Coding";
}

int BinaryTileCoding::memory()
{
    return 0;
}

MTS_NAMESPACE_END