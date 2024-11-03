#include <ds-compare/structures/binary_tile_coding.h>

MTS_NAMESPACE_BEGIN

/* ========== */
/* BinaryTile */
/* ========== */

bool BinaryTile::is_leaf() const
{
    return (this->idx_first == -1) && (this->idx_second == -1);
}

bool BinaryTile::should_split(const int depth) const
{
    if (this->data.sample_count < BinaryTileCoding::MIN_SAMPLES)
    {
        return false;
    }
    
    float diff = meandev(depth) - BinaryTileCoding::SUBDIV_THRESHOLD;
    float stderr = std::sqrt(var(depth) / this->data.sample_count);
    
    float t_own = diff / stderr;
    float t_req = TTable95::fetch(this->data.sample_count - 1);

    //std::cout << "meandev - 0.05 = diff | " << meandev(depth) << " - 0.05 = " << diff << std::endl;
    //std::cout << "stddev / sqrt(n) = stderr | " << std::sqrt(var(depth)) << " / " << std::sqrt(this->data.sample_count) << " = " << stderr << std::endl;
    //std::cout << "t_own: " << t_own << " | t_req: " << t_req << std::endl;

    return (t_own > t_req);
}

SplitDirection BinaryTile::split_direction() const
{
    return static_cast<SplitDirection>(this->data.split);
}

void BinaryTile::update_statistics(const Sample& sample)
{
    auto samples = this->data.sample_count;
    float value_mean = (samples > 0)
        ? (this->sum / samples)
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

    // Update mean deviation
    value_mean += dx / n;
    float dx2 = sample.value - value_mean;
    this->m2 += dx * dx2;
    this->diff_sum += std::abs(dx);
}

void BinaryTile::update_sum(const Sample& sample)
{
    this->sum += sample.value;
    this->data.sample_count++;
}

float BinaryTile::covar(const SplitDirection dir) const
{
    if (this->data.sample_count < 2) return 0.0f;

    float c = (dir == HORIZONTAL) ? this->cov.x : this->cov.y;
    return c / (this->data.sample_count - 1);
}

float BinaryTile::adjusted_covar(SplitDirection dir) const
{
    return std::sqrt(std::abs(covar(dir)));
}

float BinaryTile::meandev(const int depth) const
{
    if (this->data.sample_count == 0) return 0.0f;

    float area = 1.0f / (1 << depth);
    return area * this->diff_sum / this->data.sample_count;
}

float BinaryTile::var(const int depth) const
{
    if (this->data.sample_count < 2) return 0.0f;

    float area = 1.0f / (1 << depth);
    return area * area * this->m2 / (this->data.sample_count - 1);
}

float BinaryTile::mean() const
{
    if (this->data.sample_count == 0)
    {
        return 1 / (4.0f * M_PI);
    }

    return this->sum / this->data.sample_count;
}

/* ============ */
/* BinaryTiling */
/* ============ */

BinaryTile* BinaryTiling::find_tile(const Point2& uv, const Point2i& tile_dims)
{
    int unused = 0;
    return find_tile(uv, tile_dims, unused);
}

BinaryTile* BinaryTiling::find_tile(const Point2& uv, const Point2i& tile_dims, int& depth)
{
    Point2i index(
        uv.x * tile_dims.x,
        uv.y * tile_dims.y
    );

    int i = (index.y * tile_dims.x) + index.x;
    BinaryTile* curr_tile = &this->tiles.at(i);

    // Calculate boundary of current tile
    Point2 x_bounds(index.x / (Float) tile_dims.x, (index.x + 1) / (Float) tile_dims.x);
    Point2 y_bounds(index.y / (Float) tile_dims.y, (index.y + 1) / (Float) tile_dims.y);

    // Iterate through tree if necessary
    while (!curr_tile->is_leaf())
    {
        SplitDirection split_dir = curr_tile->split_direction();

        Point2& bounds = (split_dir == HORIZONTAL) ? x_bounds : y_bounds;
        Float pos = (split_dir == HORIZONTAL) ? uv.x : uv.y;
        Float split = (bounds.x + bounds.y) * 0.5;

        if (pos < split) // we're in the left or upper subtree
        {
            bounds.y = split;
            curr_tile = &this->tiles.at(curr_tile->idx_first);
        }
        else // we're in the right or lower subtree
        {
            bounds.x = split;
            curr_tile = &this->tiles.at(curr_tile->idx_second);
        }

        depth++;
    }

    return curr_tile;
}

void BinaryTiling::insert(const Sample& sample, const Point2i& tile_dims)
{
    // Find initial tile
    Point2 uv = Converter::spherical_to_uv(Point2(sample.phi, sample.theta));

    Point2 warped_pos(
        this->factors.x * (uv.x - this->start_vals.x),
        this->factors.y * (uv.y - this->start_vals.y)
    );

    int tile_depth = 0;
    BinaryTile* tile = find_tile(warped_pos, tile_dims, tile_depth);

    Sample updater;
    updater.value = sample.value;
    updater.phi = warped_pos.x;
    updater.theta = warped_pos.y;

    // Store value & update covariance
    tile->update_statistics(updater);
    tile->update_sum(updater);
    
    // Split if necessary
    if (!tile->should_split(tile_depth) || tile_depth > BinaryTileCoding::MAX_DEPTH) return;

    tile->idx_first = this->tiles.size();
    tile->idx_second = tile->idx_first + 1;

    this->tiles.push_back(BinaryTile());
    this->tiles.push_back(BinaryTile());
    
    tile->data.split = (tile->adjusted_covar(HORIZONTAL) > tile->adjusted_covar(VERTICAL))
        ? HORIZONTAL
        : VERTICAL;
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
    return;
}

RandomGen BinaryTileCoding::random = RandomGen();

Sample BinaryTileCoding::sample(Point2& pos)
{
    // [1] Pick one of the tilings with equal weight
    Float random = BinaryTileCoding::random.next1D();
    int i = random * tilings.size();
    BinaryTiling tiling = tilings.at(i);

    // [2] Use the passed-in random 2D pos to fetch the right base tile (if there are any)
    BinaryTile* curr_tile;
    Point2 x_bounds = Point2(0.0, 1.0);
    Point2 y_bounds = Point2(0.0, 1.0);

    if (tile_dims.x * tile_dims.y == 1)
    {
        curr_tile = &tiling.tiles.at(0);
    }
    else
    {
        // Calculate CDFs
        float total_mean = 0.0f;
        std::vector<float> row_means(tile_dims.y);
        for (int y = 0; y < tile_dims.y; ++y)
        {
            float row_mean = 0.0f;
            for (int x = 0; x < tile_dims.x; ++x)
            {
                int tile_i = (y * tile_dims.x) + x;
                row_mean += tiling.tiles.at(tile_i).mean();
            }
            total_mean += row_mean;
            row_means.push_back(row_mean / tile_dims.x);
        }
        total_mean /= (tile_dims.x * tile_dims.y);

        // y-dir sampling
        float sum_y = 0.0f; int y = 0;
        size_t y_len = row_means.size();
        for (y = 0; y < y_len; ++y)
        {
            sum_y += row_means.at(y) / total_mean;
            if (sum_y / y_len >= pos.y) break;
        }
        if (y == y_len) y -= 1;

        // x-dir sampling
        float sum_x = 0.0f; int x = 0;
        size_t x_len = tile_dims.x;
        for (x = 0; x < x_len; ++x)
        {
            int i = (y * x_len) + x;
            sum_x += tiling.tiles.at(i).mean() / row_means.at(y);
            if (sum_x / x_len >= pos.x) break;
        }
        if (x == x_len) x -= 1;

        curr_tile = &tiling.tiles.at((y * x_len) + x);
        x_bounds = Point2(x / (Float) x_len, (x + 1) / (Float) x_len);
        y_bounds = Point2(y / (Float) y_len, (y + 1) / (Float) y_len);
    }

    // [3] Generate random 1D sample and go deeper as long as the tile isn't a leaf
    while (!curr_tile->is_leaf())
    {
        Float random = BinaryTileCoding::random.next1D();

        SplitDirection split_dir = curr_tile->split_direction();
        Point2& bounds = (split_dir == HORIZONTAL) ? x_bounds : y_bounds;


    }

    Sample sample;
    return sample;
}

Float BinaryTileCoding::eval(Point2& pos)
{
    float total_value = 0.0f;
    for (auto& tiling : this->tilings)
    {
        Point2 warped_pos(
            tiling.factors.x * (pos.x - tiling.start_vals.x),
            tiling.factors.y * (pos.y - tiling.start_vals.y)
        );

        BinaryTile* tile = tiling.find_tile(warped_pos, this->tile_dims);
        total_value += tile->mean();
    }

    return (total_value / this->tilings.size());
}

void BinaryTileCoding::wipe()
{
    for (auto& tiling : this->tilings)
    {
        tiling = BinaryTiling();
        tiling.tiles.resize(this->tile_dims.x * this->tile_dims.y);
    }
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
    size_t size_self = sizeof(this) + sizeof(BinaryTileCoding); // base layer size
    size_t size_tilings = this->tilings.capacity() * sizeof(BinaryTiling); // #tilings * base tiling size
    
    size_t size_tiles = 0;
    for (const auto& tiling : this->tilings)
    {
        // #tiles * constant tile size
        size_tiles += tiling.tiles.capacity() * sizeof(BinaryTile);
    }

    return size_self + size_tilings + size_tiles;
}

MTS_NAMESPACE_END