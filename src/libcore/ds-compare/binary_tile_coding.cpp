#include <ds-compare/structures/binary_tile_coding.h>

MTS_NAMESPACE_BEGIN

/* ========== */
/* BinaryTile */
/* ========== */

bool BinaryTile::is_leaf() const
{
    return (this->idx_first == UINT32_MAX) && (this->idx_second == UINT32_MAX);
}

bool BinaryTile::should_split(const int depth) const
{
    if (this->sample_count < 2)
    {
        return false;
    }
    
    float diff = meandev(depth) - BinaryTileCoding::SUBDIV_THRESHOLD;
    float stderr = std::sqrt(var(depth) / this->sample_count);
    
    float t_own = diff / stderr;
    float t_req = TTable95::fetch(this->sample_count - 1);

    return (t_own > t_req);
}

SplitDirection BinaryTile::split_direction(const DepthCounter& depth) const
{
    if (this->adjusted_covar(HORIZONTAL, depth) > this->adjusted_covar(VERTICAL, depth))
    {
        return HORIZONTAL;
    }

    return VERTICAL;
}

void BinaryTile::update_statistics(const Sample& sample)
{
    auto samples = this->sample_count;
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
    this->sample_count++;
}

float BinaryTile::covar(const SplitDirection dir, const DepthCounter& counter) const
{
    if (this->sample_count < 2) return 0.0f;

    const Point2i td = BinaryTileCoding::tile_dims;
    
    float c = (dir == HORIZONTAL) 
        ? ((0.5f * this->cov.x) / (1 << counter.horizontal) / td.x)
        : (this->cov.y / (1 << counter.vertical) / td.y);

    return c / (this->sample_count - 1);
}

float BinaryTile::adjusted_covar(const SplitDirection dir, const DepthCounter& counter) const
{
    return std::sqrt(std::abs(covar(dir, counter)));
}

float BinaryTile::meandev(const int depth) const
{
    if (this->sample_count == 0) return 0.0f;

    return area(depth) * (this->diff_sum / this->sample_count);
}

float BinaryTile::var(const int depth) const
{
    if (this->sample_count < 2) return 0.0f;

    float a = area(depth);
    return a * a * (this->m2 / (this->sample_count - 1));
}

float BinaryTile::mean() const
{
    if (this->sample_count == 0)
    {
        return Epsilon;
    }

    return this->sum / this->sample_count;
}

float BinaryTile::area(const int depth) const
{
    const Point2i td = BinaryTileCoding::tile_dims;
    return 1.0f / ((1 << depth) * (td.x * td.y));
}

// Small disclaimer: I really hate that this function exists in this form and I am very sure it will be refactored in the future.
float BinaryTile::calc_visible_area_ratio(bool is_first, SplitDirection split_dir, Point2 x_bounds, Point2 y_bounds)
{
    // Return 1.0 if tile is fully within sample area
    if (!(x_bounds.x < 0 || x_bounds.y > 1 || y_bounds.x < 0 || y_bounds.y > 1))
    {
        return 1.0f;
    }

    // Return 0.0 if (sub)tile is fully outside sample area
    Point2& bounds = (split_dir == HORIZONTAL) ? x_bounds : y_bounds;
    Float halved = (bounds.x + bounds.y) * 0.5;
    if ((is_first && halved <= 0) || (!is_first && halved >= 1))
    {
        return 0.0f;
    }

    // Otherwise, calc area ratio...
    if (!is_first)
    {
        auto bound_swap = [](Point2& bounds) {
            Float temp = std::move(bounds.x);
            bounds.x = std::move(bounds.y);
            bounds.y = std::move(temp);
        };

        bound_swap(x_bounds);
        bound_swap(y_bounds);
    }

    Float x_half = (split_dir == HORIZONTAL) ? (x_bounds.x + x_bounds.y) * 0.5 : x_bounds.y;
    Float y_half = (split_dir == VERTICAL)   ? (y_bounds.x + y_bounds.y) * 0.5 : y_bounds.y;

    Point2 p1_outer(x_bounds.x, y_bounds.x);
    Point2 p2_outer(x_half, y_half);

    Point2 p1_inner = is_first
        ? Point2(std::max((Float) 0.0, x_bounds.x), std::max((Float) 0.0, y_bounds.x))
        : Point2(std::max((Float) 0.0, x_half), std::max((Float) 0.0, y_half));
    Point2 p2_inner = is_first
        ? Point2(std::min((Float) 1.0, x_half), std::min((Float) 1.0, y_half))
        : Point2(std::min((Float) 1.0, x_bounds.x), std::min((Float) 1.0, y_bounds.x));

    Float area_outer = std::abs(p2_outer.x - p1_outer.x) * std::abs(p2_outer.y - p1_outer.y);
    Float area_inner = std::abs(p2_inner.x - p1_inner.x) * std::abs(p2_inner.y - p1_inner.y);

    if (area_outer == 0) // This should never happen, but just in case...
    {
        return 0.0f;
    }

    return (area_inner / area_outer);
}

/* ============ */
/* BinaryTiling */
/* ============ */

BinaryTile& BinaryTiling::find_tile(const Point2& uv)
{
    DepthCounter unused;
    return find_tile(uv, unused);
}

BinaryTile& BinaryTiling::find_tile(const Point2& uv, DepthCounter& counter)
{
    const Point2i tile_dims = BinaryTileCoding::tile_dims;

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
        SplitDirection split_dir = curr_tile->split_direction(counter);

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

        counter.increment(split_dir);
    }

    return *curr_tile;
}

void BinaryTiling::insert(const Sample& sample)
{
    // Find initial tile
    Point2 uv = Converter::spherical_to_uv(Point2(sample.phi, sample.theta));
    Point2 warped_pos = warp_to_range(uv);

    DepthCounter counter;
    BinaryTile& tile = find_tile(warped_pos, counter);

    Sample updater;
    updater.value = sample.value;
    updater.phi = warped_pos.x;
    updater.theta = warped_pos.y;

    // Store value & update covariance
    tile.update_statistics(updater);
    tile.update_sum(updater);
    
    // Split if necessary
    if (counter.depth() > BinaryTileCoding::MAX_DEPTH || !tile.should_split(counter.depth())) return;

    tile.idx_first = this->tiles.size();
    tile.idx_second = tile.idx_first + 1;

    this->tiles.push_back(BinaryTile());
    this->tiles.push_back(BinaryTile());
}

std::pair<uint32_t, float> BinaryTiling::recurse_statistics(float& leaf_sum, BinaryTile& curr_tile, int depth)
{
    if (curr_tile.is_leaf())
    {
        leaf_sum += curr_tile.area(depth) * curr_tile.mean();
        return { curr_tile.sample_count, curr_tile.sum };
    }

    std::array<BinaryTile*, 2> children = {
        &this->tiles.at(curr_tile.idx_first),
        &this->tiles.at(curr_tile.idx_second)
    };

    for (auto child : children)
    {
        auto stats = recurse_statistics(leaf_sum, *child, depth++);

        curr_tile.sample_count += stats.first;
        curr_tile.sum += stats.second;
    }

    return { curr_tile.sample_count, curr_tile.sum };
}

Point2 BinaryTiling::warp_to_range(const Point2& uv) const
{
    return Point2(
        (uv.x - this->x_bounds.x) / (this->x_bounds.y - this->x_bounds.x),
        (uv.y - this->y_bounds.x) / (this->y_bounds.y - this->y_bounds.x)
    );
}

/* ================ */
/* BinaryTileCoding */
/* ================ */

RandomGen BinaryTileCoding::random = RandomGen();
Point2i BinaryTileCoding::tile_dims = Point2i(1, 1);
float BinaryTileCoding::SUBDIV_THRESHOLD = 0.001f;
int BinaryTileCoding::MAX_DEPTH = 10;

void BinaryTileCoding::construct(DSArguments& init_data)
{
    SAssert(init_data.btc.tilings > 0);
    SAssert(init_data.btc.tiles_x > 0 && init_data.btc.tiles_y > 0);
    SAssert(init_data.btc.max_depth > 0);

    BinaryTileCoding::tile_dims = Point2i(init_data.btc.tiles_x, init_data.btc.tiles_y);
    BinaryTileCoding::SUBDIV_THRESHOLD = init_data.btc.subdiv_threshold;
    BinaryTileCoding::MAX_DEPTH = init_data.btc.max_depth;

    this->tilings = std::vector<BinaryTiling>(init_data.btc.tilings);
    for (auto& tiling : this->tilings)
    {
        tiling.tiles = std::vector<BinaryTile>(init_data.btc.tiles_x * init_data.btc.tiles_y);
    }
}

void BinaryTileCoding::preprocess()
{
    // Calculate tiling offset
    Float tile_width = 1.0 / tile_dims.x;
    Float tile_height = 1.0 / tile_dims.y;
    auto num_tilings = this->tilings.size();

    Float overhead = 0.0;
    if (num_tilings > 1)
    {
        overhead = 1.0 / num_tilings;
    }

    Point2 offset(tile_width * overhead, tile_height * overhead);

    // Store start and end thresholds in each tiling
    for (int i = 0; i < num_tilings; ++i)
    {
        BinaryTiling& tiling = this->tilings.at(i);

        Point2 x_range(0 - (i * offset.x), 1 + ((num_tilings - 1 - i) * offset.x));
        Point2 y_range(0 - ((num_tilings - 1 - i) * offset.y), 1 + (i * offset.y));

        tiling.x_bounds = x_range;
        tiling.y_bounds = y_range;
    }
}

void BinaryTileCoding::store(std::vector<Sample>& samples)
{
    for (auto& sample : samples)
    {
        for (auto& tiling : this->tilings)
        {
            tiling.insert(sample);
        }
    }
}

void BinaryTileCoding::postprocess()
{
    for (BinaryTiling& tiling : this->tilings)
    {
        float leaf_sum = 0.0f;

        const Point2i td = BinaryTileCoding::tile_dims;
        for (int i = 0; i < td.x * td.y; ++i)
        {
            BinaryTile& tile = tiling.tiles.at(i);
            tiling.recurse_statistics(leaf_sum, tile, 0);
        }

        tiling.leaf_sum = leaf_sum;
    }
}

Sample BinaryTileCoding::sample(Point2& pos)
{
    // [1] Pick one of the tilings with equal weight
    Float random = BinaryTileCoding::random.next1D();
    
    int i = random * this->tilings.size();
    /*Float total_leaf_sum = 0;
    for (const auto& tiling : this->tilings)
    {
        total_leaf_sum += tiling.leaf_sum;
    }

    Float sum_leafs = 0; int i = 0;
    for (i = 0; i < this->tilings.size(); ++i)
    {
        sum_leafs += this->tilings.at(i).leaf_sum / total_leaf_sum;
        if (sum_leafs >= random) break;
    }*/
    BinaryTiling& tiling = tilings.at(i);

    // [2] Use the passed-in random 2D pos to fetch the right base tile (if there are any)
    BinaryTile* curr_tile;
    Point2 x_bounds = tiling.x_bounds;
    Point2 y_bounds = tiling.y_bounds;
    
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
            row_means.at(y) = (row_mean / tile_dims.x);
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
        Float dx = tiling.x_bounds.y - tiling.x_bounds.x;
        Float dy = tiling.y_bounds.y - tiling.y_bounds.x;
        x_bounds = Point2(tiling.x_bounds.x + (x / (Float) x_len) * dx, tiling.x_bounds.x + ((x + 1) / (Float) x_len) * dx);
        y_bounds = Point2(tiling.y_bounds.x + (y / (Float) y_len) * dy, tiling.y_bounds.x + ((y + 1) / (Float) y_len) * dy);
    }

    // [3] Generate random 1D sample and go deeper as long as the tile isn't a leaf
    DepthCounter counter;
    Float area_mult = 1.0;
    while (!curr_tile->is_leaf())
    {
        Float random = BinaryTileCoding::random.next1D();

        SplitDirection split_dir = curr_tile->split_direction(counter);
        Point2& bounds = (split_dir == HORIZONTAL) ? x_bounds : y_bounds;
        Float halved = (bounds.x + bounds.y) * 0.5;

        BinaryTile* first = &tiling.tiles.at(curr_tile->idx_first);
        BinaryTile* second = &tiling.tiles.at(curr_tile->idx_second);

        Float first_area_mult = BinaryTile::calc_visible_area_ratio(true, split_dir, x_bounds, y_bounds);
        Float second_area_mult = BinaryTile::calc_visible_area_ratio(false, split_dir, x_bounds, y_bounds);
        Float first_mean = first_area_mult * first->mean();
        Float second_mean = second_area_mult * second->mean();

        Float split = first_mean / (first_mean + second_mean);

        if (random < split)
        {
            bounds.y = halved;
            curr_tile = first;
            area_mult = first_area_mult;
        }
        else
        {
            bounds.x = halved;
            curr_tile = second;
            area_mult = second_area_mult;
        }
        
        counter.increment(split_dir);
    }

    Point2 rng = BinaryTileCoding::random.next2D();
    Point2 coords(
        x_bounds.x + rng.x * (x_bounds.y - x_bounds.x),
        y_bounds.x + rng.y * (y_bounds.y - y_bounds.x)
    );

    if (x_bounds.x < 0)
        coords.x = ((coords.x - x_bounds.x) * x_bounds.y) / (x_bounds.y - x_bounds.x);

    if (y_bounds.x < 0)
        coords.y = ((coords.y - y_bounds.x) * y_bounds.y) / (y_bounds.y - y_bounds.x);

    if (x_bounds.y > 1)
        coords.x = x_bounds.x + ((coords.x - x_bounds.x) * (1.0 - x_bounds.x)) / (x_bounds.y - x_bounds.x);

    if (y_bounds.y > 1)
        coords.y = y_bounds.x + ((coords.y - y_bounds.x) * (1.0 - y_bounds.x)) / (y_bounds.y - y_bounds.x);

    float area = curr_tile->area(counter.depth());
    float prob = area_mult * area * (curr_tile->mean() / tiling.leaf_sum);

    Sample sample = {
        .value = 0,
        .pdf = prob,
        .theta = (coords.y * M_PI),
        .phi = (coords.x * 2 * M_PI)
    };
    return sample;
}

Float BinaryTileCoding::eval(Point2& pos)
{
    float total_value = 0.0f;
    for (auto& tiling : this->tilings)
    {
        Point2 warped_pos = tiling.warp_to_range(pos);
        BinaryTile& tile = tiling.find_tile(warped_pos);
        total_value += tile.mean();
    }

    return (total_value / this->tilings.size());
}

void BinaryTileCoding::wipe()
{
    for (auto& tiling : this->tilings)
    {
        tiling = BinaryTiling();
        tiling.tiles = std::vector<BinaryTile>(tile_dims.x * tile_dims.y);
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