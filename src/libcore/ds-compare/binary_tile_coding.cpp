#include <ds-compare/structures/binary_tile_coding.h>

MTS_NAMESPACE_BEGIN

/* ========== */
/* BinaryTile */
/* ========== */

bool BinaryTile::is_leaf() const
{
    return (this->idx_first == UINT32_MAX) && (this->idx_second == UINT32_MAX);
}

bool BinaryTile::should_split(const TileTracker& tracker) const
{
    if (this->sample_count < 2)
    {
        return false;
    }
    
    float diff = meandev(tracker) - BinaryTileCoding::SUBDIV_THRESHOLD;
    float stderr = std::sqrt(var(tracker) / this->sample_count);
    
    float t_own = diff / stderr;
    float t_req = TTable95::fetch(this->sample_count - 1);

    return (t_own > t_req);
}

Direction BinaryTile::split_direction(const TileTracker& tracker) const
{
    if (this->adjusted_covar(HORIZONTAL, tracker) > this->adjusted_covar(VERTICAL, tracker))
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

float BinaryTile::covar(const Direction dir, const TileTracker& tracker) const
{
    if (this->sample_count < 2) return 0.0f;

    auto x_norm = tracker.clamped(HORIZONTAL);
    if (x_norm.x == x_norm.y) return 0.0f;

    auto y_norm = tracker.clamped(VERTICAL);
    if (y_norm.x == y_norm.y) return 0.0f;
    
    float c = (dir == HORIZONTAL) ? this->cov.x : this->cov.y;

    return c / (this->sample_count - 1);
}

float BinaryTile::adjusted_covar(const Direction dir, const TileTracker& tracker) const
{
    return std::sqrt(std::abs(covar(dir, tracker)));
}

float BinaryTile::meandev(const TileTracker& tracker) const
{
    if (this->sample_count == 0) return 0.0f;

    return area(tracker) * (this->diff_sum / this->sample_count);
}

float BinaryTile::var(const TileTracker& tracker) const
{
    if (this->sample_count < 2) return 0.0f;

    float a = area(tracker);
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

float BinaryTile::area(const TileTracker& tracker) const
{
    auto x_clamped = tracker.clamped(HORIZONTAL);
    if (x_clamped.x == x_clamped.y) return 0.0f;

    auto y_clamped = tracker.clamped(VERTICAL);
    if (y_clamped.x == y_clamped.y) return 0.0f;
    
    Float d_theta = std::cos(M_PI * y_clamped.x) - std::cos(M_PI * y_clamped.y);
    Float d_phi = 2 * M_PI * (x_clamped.y - x_clamped.x);

    return d_theta * d_phi;
}

float BinaryTile::visible_area_perc(bool before_split, Direction split_dir, Point2 x_bounds, Point2 y_bounds)
{
    // Return 1.0 if tile is fully within sample area
    if (!(x_bounds.x < 0 || x_bounds.y > 1 || y_bounds.x < 0 || y_bounds.y > 1))
    {
        return 1.0f;
    }

    // Return 0.0 if tile is fully outside sample area
    Point2& bounds = (split_dir == HORIZONTAL) ? x_bounds : y_bounds;
    Float halved = (bounds.x + bounds.y) * 0.5;
    if ((before_split && halved <= 0) || (!before_split && halved >= 1))
    {
        return 0.0f;
    }

    // Otherwise, calc area ratio...
    if (!before_split)
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

    Point2 p1_inner = before_split
        ? Point2(std::max((Float) 0.0, x_bounds.x), std::max((Float) 0.0, y_bounds.x))
        : Point2(std::max((Float) 0.0, x_half), std::max((Float) 0.0, y_half));
    Point2 p2_inner = before_split
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

BinaryTile& BinaryTiling::find_tile(const Point2& pos)
{
    TileTracker unused;
    return find_tile(pos, unused);
}

BinaryTile& BinaryTiling::find_tile(const Point2& pos, TileTracker& tracker)
{
    const Point2i td = BinaryTileCoding::tile_dims;

    // Find right base tile by mapping the pos within the range of the tiling to the range [0, 1]
    auto warped_pos = warp_to_range(pos);
    Point2i bt_pos(
        warped_pos.x * td.x,
        warped_pos.y * td.y
    );

    int i = (bt_pos.y * td.x) + bt_pos.x;
    BinaryTile* curr_tile = &this->tiles.at(i);

    // Calculate boundary of current tile
    float x_step = (this->x_bounds.y - this->x_bounds.x) / td.x;
    float y_step = (this->y_bounds.y - this->y_bounds.x) / td.y;

    Point2 x_bounds(
        this->x_bounds.x + (bt_pos.x * x_step),
        this->x_bounds.x + ((bt_pos.x + 1) * x_step)
    );
    Point2 y_bounds(
        this->y_bounds.x + (bt_pos.y * y_step),
        this->y_bounds.x + (bt_pos.y + 1) * y_step
    );

    // Iterate through tree if necessary
    while (!curr_tile->is_leaf())
    {
        Direction split_dir = curr_tile->split_direction(tracker);

        Point2& bounds = (split_dir == HORIZONTAL) ? x_bounds : y_bounds;
        Float local = (split_dir == HORIZONTAL) ? pos.x : pos.y;
        Float split = (bounds.x + bounds.y) * 0.5;

        if (local < split) // we're in the left or upper subtree
        {
            bounds.y = split;
            curr_tile = &this->tiles.at(curr_tile->idx_first);
        }
        else // we're in the right or lower subtree
        {
            bounds.x = split;
            curr_tile = &this->tiles.at(curr_tile->idx_second);
        }

        tracker.increment(split_dir);
        tracker.side(local < split);
    }

    tracker.boundaries(x_bounds, y_bounds);
    return *curr_tile;
}

void BinaryTiling::insert(const Sample& sample)
{
    // Find initial tile
    TileTracker tracker;
    Point2 uv = Converter::spherical_to_uv(Point2(sample.phi, sample.theta));
    BinaryTile& tile = find_tile(uv, tracker);

    Sample updater;
    updater.value = sample.value;
    updater.phi = uv.x;
    updater.theta = uv.y;

    // Store value & update covariance
    tile.update_statistics(sample);
    tile.update_sum(sample);
    
    // Split if necessary
    if (tracker.depth() > BinaryTileCoding::MAX_DEPTH || !tile.should_split(tracker)) return;

    tile.idx_first = this->tiles.size();
    tile.idx_second = tile.idx_first + 1;

    this->tiles.push_back(BinaryTile());
    this->tiles.push_back(BinaryTile());
}

std::pair<uint32_t, float> BinaryTiling::recurse_statistics(BinaryTile& curr_tile, TileTracker tracker, float& leaf_sum)
{
    if (curr_tile.is_leaf())
    {
        float area = curr_tile.area(tracker);
        float mu = curr_tile.mean();
        leaf_sum += area * mu;

        return { curr_tile.sample_count, curr_tile.sum };
    }

    std::array<BinaryTile*, 2> children = {
        &this->tiles.at(curr_tile.idx_first),
        &this->tiles.at(curr_tile.idx_second)
    };

    Direction split_dir = curr_tile.split_direction(tracker);
    tracker.increment(split_dir);

    for (int i = 0; i < 2; ++i)
    {
        TileTracker c_tracker = tracker;
        c_tracker.side(i == 0);

        Point2& bounds = (c_tracker.last == HORIZONTAL) ? c_tracker.x_bounds : c_tracker.y_bounds;
        Float& value = c_tracker.before_split ? bounds.y : bounds.x;
        value = (bounds.x + bounds.y) * 0.5;

        BinaryTile* child = children.at(i);
        auto stats = recurse_statistics(*child, c_tracker, leaf_sum);

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

Point2i BinaryTiling::base_tile_pos_from_cdf(const Point2& pos) const
{
    Point2i td = BinaryTileCoding::tile_dims;

    // Calculate CDFs
    float total_mean = 0.0f;
    std::vector<float> row_means(td.y);
    for (int y = 0; y < td.y; ++y)
    {
        float row_mean = 0.0f;
        for (int x = 0; x < td.x; ++x)
        {
            int tile_i = (y * td.x) + x;
            row_mean += this->tiles.at(tile_i).mean();
        }
        total_mean += row_mean;
        row_means.at(y) = (row_mean / td.x);
    }
    total_mean /= (td.x * td.y);

    // y-dir sampling
    float sum_y = 0.0f; int y = 0;
    size_t y_len = row_means.size();
    for (y; y < y_len; ++y)
    {
        sum_y += row_means.at(y) / total_mean;
        if (sum_y / y_len >= pos.y) break;
    }
    if (y == y_len) y -= 1;

    // x-dir sampling
    float sum_x = 0.0f; int x = 0;
    size_t x_len = td.x;
    for (x; x < x_len; ++x)
    {
        int i = (y * x_len) + x;
        sum_x += this->tiles.at(i).mean() / row_means.at(y);
        if (sum_x / x_len >= pos.x) break;
    }
    if (x == x_len) x -= 1;

    return Point2i(x, y);
}

float BinaryTiling::pdf(const Point2& pos)
{
    TileTracker tracker;
    BinaryTile& tile = find_tile(pos, tracker);

    float area = tile.area(tracker);
    float mu = tile.mean();

    return area * mu;
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
    auto num_tilings = this->tilings.size();

    Float overhead = 0.0;
    if (num_tilings > 1)
    {
        overhead = 1.0 / num_tilings;
    }

    // Store start and end thresholds in each tiling
    for (int i = 0; i < num_tilings; ++i)
    {
        BinaryTiling& tiling = this->tilings.at(i);

        Point2 x_range(0 - (i * overhead), 1 + ((num_tilings - 1 - i) * overhead));
        Point2 y_range(0 - ((num_tilings - 1 - i) * overhead), 1 + (i * overhead));

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
        const Point2i td = BinaryTileCoding::tile_dims;
        Point2 step(
            (tiling.x_bounds.y - tiling.x_bounds.x) / td.x,
            (tiling.y_bounds.y - tiling.y_bounds.x) / td.y
        );

        for (int i = 0; i < td.x * td.y; ++i)
        {
            TileTracker tracker;

            Point2i pos(i % td.x, i / td.x);
            Point2 local_x_bounds(tiling.x_bounds.x + step.x * pos.x, tiling.x_bounds.x + step.x * (pos.x + 1));
            Point2 local_y_bounds(tiling.y_bounds.x + step.y * pos.y, tiling.y_bounds.x + step.y * (pos.y + 1));
            tracker.boundaries(local_x_bounds, local_y_bounds);
            
            BinaryTile& tile = tiling.tiles.at(i);
            tiling.recurse_statistics(tile, tracker, this->leaf_sum);
        }
    }
}

Sample BinaryTileCoding::sample(Point2& pos)
{
    Sample sample = {
           .value = 0,
           .pdf = 0.01,
           .theta = 0,
           .phi = 0
    };
    return sample;

    //int tiling_count = this->tilings.size();

    //// [1] Pick one of the tilings with equal weight
    //Float random = BinaryTileCoding::random.next1D();
    //int i = random * tiling_count;
    //BinaryTiling& tiling = this->tilings.at(i);

    //// [2] If there is more than one base tile, use the passed-in (ideally random) 2D position to fetch the right one
    //BinaryTile* curr_tile = &tiling.tiles.at(0);
    //Point2 x_bounds = tiling.x_bounds;
    //Point2 y_bounds = tiling.y_bounds;
    //
    //if (tile_dims.x * tile_dims.y > 1)
    //{
    //    auto bt_pos = tiling.base_tile_pos_from_cdf(pos);
    //    curr_tile = &tiling.tiles.at((bt_pos.y * tile_dims.x) + bt_pos.x);

    //    Float x_step = (tiling.x_bounds.y - tiling.x_bounds.x) / tile_dims.x;
    //    Float y_step = (tiling.y_bounds.y - tiling.y_bounds.x) / tile_dims.y;

    //    x_bounds = Point2(
    //        tiling.x_bounds.x + (bt_pos.x * x_step),
    //        tiling.x_bounds.x + ((bt_pos.x + 1) * x_step)
    //    );
    //    y_bounds = Point2(
    //        tiling.y_bounds.x + (bt_pos.y * y_step),
    //        tiling.y_bounds.x + ((bt_pos.y + 1) * y_step)
    //    );
    //}

    //// [3] Generate random 1D sample and go deeper as long as the tile isn't a leaf
    //TileTracker tracker;
    //Float visible_perc = 1.0;
    //while (!curr_tile->is_leaf())
    //{
    //    Float random = BinaryTileCoding::random.next1D();

    //    Direction split_dir = curr_tile->split_direction(tracker);
    //    Point2& bounds = (split_dir == HORIZONTAL) ? x_bounds : y_bounds;
    //    Float halved = (bounds.x + bounds.y) * 0.5;

    //    BinaryTile* first = &tiling.tiles.at(curr_tile->idx_first);
    //    BinaryTile* second = &tiling.tiles.at(curr_tile->idx_second);

    //    Float area_first = BinaryTile::visible_area_perc(true, split_dir, x_bounds, y_bounds);
    //    Float area_second = BinaryTile::visible_area_perc(false, split_dir, x_bounds, y_bounds);
    //    Float first_mean = area_first * first->mean();
    //    Float second_mean = area_second * second->mean();

    //    Float split = first_mean / (first_mean + second_mean);

    //    if (random < split)
    //    {
    //        bounds.y = halved;
    //        curr_tile = first;
    //        visible_perc = area_first;
    //    }
    //    else
    //    {
    //        bounds.x = halved;
    //        curr_tile = second;
    //        visible_perc = area_second;
    //    }

    //    tracker.increment(split_dir);
    //}

    //// [4] Generate random sample once in a leaf tile
    //Point2 rng = BinaryTileCoding::random.next2D();
    //Point2 coords(
    //    x_bounds.x + rng.x * (x_bounds.y - x_bounds.x),
    //    y_bounds.x + rng.y * (y_bounds.y - y_bounds.x)
    //);

    //if (x_bounds.x < 0)
    //    coords.x = ((coords.x - x_bounds.x) * x_bounds.y) / (x_bounds.y - x_bounds.x);

    //if (y_bounds.x < 0)
    //    coords.y = ((coords.y - y_bounds.x) * y_bounds.y) / (y_bounds.y - y_bounds.x);

    //if (x_bounds.y > 1)
    //    coords.x = x_bounds.x + ((coords.x - x_bounds.x) * (1.0 - x_bounds.x)) / (x_bounds.y - x_bounds.x);

    //if (y_bounds.y > 1)
    //    coords.y = y_bounds.x + ((coords.y - y_bounds.x) * (1.0 - y_bounds.x)) / (y_bounds.y - y_bounds.x);

    //// [5] Calculate PDF value based on position
    //float area = curr_tile->area(tracker.depth());
    //float mu = curr_tile->mean();
    //float prob = visible_perc * area * mu;

    //for (int ti = 0; ti < tiling_count; ++ti)
    //{
    //    if (ti == i) continue;
    //    prob += this->tilings.at(ti).pdf(coords);
    //}

    //prob /= this->leaf_sum;

    //Sample sample = {
    //    .value = 0,
    //    .pdf = prob,
    //    .theta = (coords.y * M_PI),
    //    .phi = (coords.x * 2 * M_PI)
    //};
    //return sample;
}

Float BinaryTileCoding::eval(Point2& pos)
{
    float prob = 0.0f;
    for (auto& tiling : this->tilings)
    {
        prob += tiling.pdf(pos);
    }
    
    return (prob / this->leaf_sum);
}

void BinaryTileCoding::wipe()
{
    this->leaf_sum = 0.0f;
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