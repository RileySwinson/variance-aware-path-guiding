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
    if (this->sample_count < BinaryTileCoding::MIN_SAMPLES)
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

/*BinaryTile BinaryTile::stat_copy_of(const BinaryTile& other)
{
    BinaryTile tile;
    tile.cov = other.cov;
    tile.sample_mean = other.sample_mean;
    tile.m2 = other.m2;
    return tile;
}*/

float BinaryTile::covar(const SplitDirection dir, const DepthCounter& counter) const
{
    if (this->sample_count < 2) return 0.0f;
    
    float c = (dir == HORIZONTAL) 
        ? ((0.5f * this->cov.x) / (1 << counter.horizontal))
        : (this->cov.y / (1 << counter.vertical));

    return c / (this->sample_count - 1);
}

float BinaryTile::adjusted_covar(const SplitDirection dir, const DepthCounter& counter) const
{
    return std::sqrt(std::abs(covar(dir, counter)));
}

float BinaryTile::meandev(const int depth) const
{
    if (this->sample_count == 0) return 0.0f;

    float area = 1.0f / (1 << depth);
    return area * (this->diff_sum / this->sample_count);
}

float BinaryTile::var(const int depth) const
{
    if (this->sample_count < 2) return 0.0f;

    float area = 1.0f / (1 << depth);
    return area * area * (this->m2 / (this->sample_count - 1));
}

float BinaryTile::mean() const
{
    if (this->sample_count == 0)
    {
        return Epsilon;
    }

    return this->sum / this->sample_count;
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

float BinaryTiling::calc_leaf_sum()
{
    typedef std::pair<BinaryTile*, int> d_tile;
    auto create_tile = [this](int idx, int depth) { return d_tile(&this->tiles.at(idx), depth); };

    std::stack<d_tile> tile_storage;
    tile_storage.push(create_tile(0, 0));

    float leaf_sum = 0.0f;
    while (!tile_storage.empty())
    {
        d_tile& dt = tile_storage.top();
        tile_storage.pop();

        BinaryTile* curr_tile = dt.first;
        int depth = dt.second;

        if (curr_tile->is_leaf())
        {
            float area = 1.0f / (1 << depth);
            leaf_sum += curr_tile->mean() * area;
            continue;
        }

        tile_storage.push(create_tile(curr_tile->idx_first, depth++));
        tile_storage.push(create_tile(curr_tile->idx_second, depth++));
    }

    return leaf_sum;
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
int BinaryTileCoding::MIN_SAMPLES = 1000;
int BinaryTileCoding::MAX_DEPTH = 10;

void BinaryTileCoding::construct(DSArguments& init_data)
{
    SAssert(init_data.btc.tilings > 0);
    SAssert(init_data.btc.tiles_x > 0 && init_data.btc.tiles_y > 0);
    SAssert(init_data.btc.max_depth > 0);

    BinaryTileCoding::tile_dims = Point2i(init_data.btc.tiles_x, init_data.btc.tiles_y);
    BinaryTileCoding::SUBDIV_THRESHOLD = init_data.btc.subdiv_threshold;
    BinaryTileCoding::MIN_SAMPLES = init_data.btc.min_tile_samples;
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
        tiling.area_leaf_sum = tiling.calc_leaf_sum();
    }
}

Sample BinaryTileCoding::sample(Point2& pos)
{
    // [1] Pick one of the tilings with equal weight
    Float random = BinaryTileCoding::random.next1D();
    int i = random * tilings.size();
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
    DepthCounter counter;
    while (!curr_tile->is_leaf())
    {
        Float random = BinaryTileCoding::random.next1D();

        SplitDirection split_dir = curr_tile->split_direction(counter);
        Point2& bounds = (split_dir == HORIZONTAL) ? x_bounds : y_bounds;

        BinaryTile* first = &tiling.tiles.at(curr_tile->idx_first);
        BinaryTile* second = &tiling.tiles.at(curr_tile->idx_second);

        Float split = first->mean() / (first->mean() + second->mean());
        Float halved = (bounds.x + bounds.y) * 0.5;

        if (random < split)
        {
            bounds.y = halved;
            curr_tile = first;
        }
        else
        {
            bounds.x = halved;
            curr_tile = second;
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

    float area = 1.0f / (1 << counter.depth());
    float prob = area * (curr_tile->mean() / tiling.area_leaf_sum);

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