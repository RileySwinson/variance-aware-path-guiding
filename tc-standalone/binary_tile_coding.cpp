#include "binary_tile_coding.h"

TC_NAMESPACE_BEGIN

/* ========== */
/* BinaryTile */
/* ========== */

void BinaryTile::update(const Sample& sample)
{
    if (!is_leaf()) return;

    auto samples = this->data.sample_count;
    float value_mean = (samples > 0) ? (this->sum / samples) : 0;
    const float contribution = sample.value / sample.p;
    auto n = samples + 1;

    float dx = contribution - value_mean;
    auto dy_x = [&]() { return sample.phi - this->leaf.sample_mean.x; };
    auto dy_y = [&]() { return sample.theta - this->leaf.sample_mean.y; };

    // Update x covariance
    this->leaf.sample_mean.x += dy_x() / n;
    this->leaf.cov.x += dx * dy_x() * std::sin(sample.theta);

    // Update y covariance
    this->leaf.sample_mean.y += dy_y() / n;
    this->leaf.cov.y += dx * dy_y();

    // Update mean deviation
    value_mean += dx / n;
    float dx2 = contribution - value_mean;
    float diff_squared = dx * dx2;
    this->leaf.m2 += diff_squared;
    this->leaf.ad_sum += std::sqrt(diff_squared);

    // Update sample count & sum
    this->data.sample_count++;
    this->sum += contribution;
}

Direction BinaryTile::split_direction(const BTTracker& tracker) const
{
    auto x_norm = tracker.clamped(Horizontal);
    auto y_norm = tracker.clamped(Vertical);

    const float covar_x = (x_norm.x != x_norm.y) ? std::abs(this->leaf.cov.x) : 0.0f;
    const float covar_y = (y_norm.x != y_norm.y) ? std::abs(this->leaf.cov.y) : 0.0f;

    if (covar_x > covar_y)
    {
        return Horizontal;
    }

    return Vertical;
}

bool BinaryTile::should_split() const
{
    auto sample_count = this->data.sample_count;
    if (sample_count < 2)
    {
        return false;
    }

    // We assume that the population mean is 0.
    float std_err = std::sqrt(var() / sample_count);

    float t_own = meandev() / std_err;
    float t_req = TTable::fetch(BinaryTileCoding::CI, sample_count - 1);

    return (t_own > t_req);
}

float BinaryTile::meandev() const
{
    auto sample_count = this->data.sample_count;
    if (sample_count == 0)
    {
        return 0.0f;
    }

    return (this->leaf.ad_sum / sample_count);
}

float BinaryTile::var() const
{
    auto sample_count = this->data.sample_count;
    if (sample_count < 2)
    {
        return 0.0f;
    }

    return (this->leaf.m2 / (sample_count - 1));
}

float BinaryTile::area(const BTTracker& tracker, bool transform)
{
    auto x_clamped = tracker.clamped(Horizontal);
    if (x_clamped.x == x_clamped.y) return 0.0f;

    auto y_clamped = tracker.clamped(Vertical);
    if (y_clamped.x == y_clamped.y) return 0.0f;

    if (transform)
    {
        y_clamped.x = BinaryTiling::transform(y_clamped.x);
        y_clamped.y = BinaryTiling::transform(y_clamped.y);
    }

    float d_theta = std::cos(M_PI * y_clamped.x) - std::cos(M_PI * y_clamped.y);
    float d_phi = 2 * M_PI * (x_clamped.y - x_clamped.x);

    return d_theta * d_phi;
}

bool BinaryTile::is_leaf() const
{
    return this->data.tile_type == TileType::Leaf;
}

/* ============ */
/* BinaryTiling */
/* ============ */

BinaryTile& BinaryTiling::find_tile(const Pair<float>& pos, bool find_empty)
{
    BTTracker unused;
    return find_tile(pos, unused, find_empty);
}

BinaryTile& BinaryTiling::find_tile(const Pair<float>& pos, BTTracker& tracker, bool find_empty)
{
    const Pair<int> td = BinaryTileCoding::TILE_DIMS;

    // Find right base tile by mapping the pos within the range of the tiling to the range [0, 1]
    Pair<float> pos_warped(pos.x, transform(pos.y, true));
    auto x_warped = Util::map(pos_warped.x).from({ this->x_bounds.x, this->x_bounds.y }).to({ 0.0, 1.0 });
    auto y_warped = Util::map(pos_warped.y).from({ this->y_bounds.x, this->y_bounds.y }).to({ 0.0, 1.0 });

    Pair<int> bt_pos(
        (int)(x_warped * td.x),
        (int)(y_warped * td.y)
    );

    int i = (bt_pos.y * td.x) + bt_pos.x;
    BinaryTile* curr_tile = &this->tiles[i];

    // Calculate planar boundary of current tile
    auto bounds = base_tile_bounds(bt_pos.x, bt_pos.y);
    Pair<float> x_bounds = bounds.first;
    Pair<float> y_bounds = bounds.second;

    // Iterate through tree if necessary
    while (!curr_tile->is_leaf())
    {
        InternalNode& parent = curr_tile->internal;
        Direction split_dir = parent.split_direction;
        bool h_split = (split_dir == Horizontal);

        Pair<float>& bounds = h_split ? x_bounds : y_bounds;
        float local = h_split ? pos.x : pos.y;
        float split = (bounds.x + bounds.y) * 0.5f;

        if (!find_empty && children_empty(*curr_tile))
        {
            break;
        }

        float split_pos = h_split ? split : transform(split);
        bool first_child = local < split_pos;

        if (first_child)
        {
            bounds.y = split; // we're in the left or upper subtree
        }
        else
        {
            bounds.x = split; // we're in the right or lower subtree
        }

        tracker.increment();
        tracker.sum = parent.power[!first_child];

        auto index = parent.children[!first_child];
        curr_tile = &this->tiles[index];
    }

    tracker.set_boundaries(x_bounds, y_bounds);
    return *curr_tile;
}

void BinaryTiling::insert(const Sample& sample, uint32_t& split_count)
{
    // Find initial tile
    BTTracker tracker;
    Pair<float> uv = Util::spherical_to_uv(Pair<float>(sample.phi, sample.theta));
    BinaryTile& tile = find_tile(uv, tracker);

    // Store value & update covariance
    tile.update(sample);

    // Only split if necessary
    bool below_max_depth = (tracker.depth < BinaryTileCoding::MAX_DEPTH);
    bool below_max_splits = (split_count < BinaryTileCoding::MAX_SPLITS);

    if (below_max_depth && below_max_splits && tile.should_split())
    {
        uint32_t tiles = this->tiles.size();

        InternalNode parent;
        parent.children = { tiles, tiles + 1 };
        parent.split_direction = tile.split_direction(tracker);

        tile.internal = parent;
        tile.data.tile_type = TileType::Internal;

        BinaryTile child(Leaf);
        this->tiles.push_back(child);
        this->tiles.push_back(child);

        split_count++;
    }
}

float BinaryTiling::recurse_statistics(BinaryTile& curr_tile, BTTracker tracker, float& leaf_sum)
{
    if (curr_tile.is_leaf() || children_empty(curr_tile))
    {
        float power = curr_tile.sum * curr_tile.area(tracker, true);
        leaf_sum += power;

        return power;
    }

    InternalNode& parent = curr_tile.internal;
    parent.power.fill(0.0f);

    for (int i = 0; i < 2; ++i)
    {
        BinaryTile& child = this->tiles[parent.children[i]];
        if (child.sum == 0) continue;

        BTTracker c_tracker = tracker;
        c_tracker.increment();

        Pair<float>& bounds = (parent.split_direction == Horizontal)
            ? c_tracker.x_bounds
            : c_tracker.y_bounds;

        float& value = (i == 0) ? bounds.y : bounds.x;
        value = (bounds.x + bounds.y) * 0.5;

        parent.power[i] += recurse_statistics(child, c_tracker, leaf_sum);
    }

    return parent.power[0] + parent.power[1];
}

std::pair<Pair<float>, Pair<float>> BinaryTiling::base_tile_bounds(int x, int y) const
{
    Pair<int> td = BinaryTileCoding::TILE_DIMS;

    float x_step = (this->x_bounds.y - this->x_bounds.x) / td.x;
    float y_step = (this->y_bounds.y - this->y_bounds.x) / td.y;

    Pair<float> x_bounds(
        this->x_bounds.x + (x * x_step),
        this->x_bounds.x + ((x + 1) * x_step)
    );
    Pair<float> y_bounds(
        this->y_bounds.x + (y * y_step),
        this->y_bounds.x + ((y + 1) * y_step)
    );

    return { x_bounds, y_bounds };
}

Pair<int> BinaryTiling::base_tile_pos(const Pair<float>& pos)
{
    // Marginal sampling (y-dir)
    size_t y = 0;
    float sum_y = 0.0f;

    size_t y_len = BinaryTileCoding::TILE_DIMS.y;
    for (; y < y_len; ++y)
    {
        sum_y += this->row_means[y];

        if (sum_y / (y_len * this->total_power) >= pos.y)
        {
            break;
        }
    }

    if (y == y_len) y -= 1;

    // Conditional sampling (x-dir)
    size_t x = 0;
    float sum_x = 0.0f;

    size_t x_len = BinaryTileCoding::TILE_DIMS.x;
    for (; x < x_len; ++x)
    {
        int i = (y * x_len) + x;
        auto& curr = this->tiles[i].internal;

        float power = curr.power[0] + curr.power[1];
        sum_x += power * this->tile_areas[i];

        if (sum_x / (x_len * this->row_means[y]) >= pos.x)
        {
            break;
        }
    }

    if (x == x_len) x -= 1;

    return Pair<int>(x, y);
}

inline bool BinaryTiling::children_empty(BinaryTile& tile) const
{
    if (tile.is_leaf())
    {
        return false;
    }

    return (this->tiles[tile.internal.children[0]].sum == 0) && (this->tiles[tile.internal.children[1]].sum == 0);
}

template <typename T>
T BinaryTiling::transform(T value, bool inverse, const std::function<T(T)>& f)
{
    if (f) return f(value);
    return Transformation::f<T>(BinaryTileCoding::TRANSFORM_MODE, inverse)(value);
}

/* ================ */
/* BinaryTileCoding */
/* ================ */

RandomGen BinaryTileCoding::random = RandomGen();
Pair<int> BinaryTileCoding::TILE_DIMS = Pair<int>(1, 1);
float BinaryTileCoding::TILING_EXCESS = 0.1f;
uint32_t BinaryTileCoding::MAX_DEPTH = 10;
uint32_t BinaryTileCoding::MAX_SPLITS = UINT32_MAX;
TTable::CI BinaryTileCoding::CI = TTable::CI::P950;
Transformation::Domain BinaryTileCoding::TRANSFORM_MODE = Transformation::Domain::Spherical;

void BinaryTileCoding::construct(const BTCArguments& args)
{
    assert(args.tilings > 0);
    assert(args.tiles.x > 0 && args.tiles.y > 0);
    assert(args.max_depth >= 0);
    assert(args.eagerness > -1 && args.eagerness < 5);

    BinaryTileCoding::TILE_DIMS = args.tiles;
    BinaryTileCoding::TILING_EXCESS = args.excess;
    BinaryTileCoding::MAX_DEPTH = args.max_depth;
    BinaryTileCoding::MAX_SPLITS = args.max_splits;
    BinaryTileCoding::CI = static_cast<TTable::CI>(args.eagerness);
    BinaryTileCoding::TRANSFORM_MODE = args.transformation;

    this->tilings = std::vector<BinaryTiling>(args.tilings);
}

void BinaryTileCoding::preprocess()
{
    // Allocate tiles space
    int base_tiles = TILE_DIMS.x * TILE_DIMS.y;
    int max_cap = (1 << MAX_DEPTH) - 1;

    for (auto& tiling : this->tilings)
    {
        tiling.row_means.resize(TILE_DIMS.y);
        tiling.tile_areas.resize(base_tiles);

        tiling.tiles.reserve(base_tiles * max_cap);

        for (int i = 0; i < base_tiles; ++i)
        {
            tiling.tiles.push_back(BinaryTile(Leaf));
        }
    }

    // Calculate tiling offset
    auto tiling_count = this->tilings.size();

    float shift = 0.0;
    if (tiling_count > 1)
    {
        shift = TILING_EXCESS / (tiling_count - 1);
    }

    // Store start and end thresholds in each tiling
    for (size_t i = 0; i < tiling_count; ++i)
    {
        BinaryTiling& tiling = this->tilings[i];

        Pair<float> x_range(0 - (i * shift), 1 + ((tiling_count - 1 - i) * shift));
        Pair<float> y_range(0 - ((tiling_count - 1 - i) * shift), 1 + (i * shift));

        tiling.x_bounds = x_range;
        tiling.y_bounds = y_range;
    }
}

void BinaryTileCoding::store(Sample& sample)
{
    for (auto& tiling : this->tilings)
    {
        tiling.insert(sample, this->split_count);
    }
}

void BinaryTileCoding::build()
{
    const Pair<int> td = BinaryTileCoding::TILE_DIMS;
    this->leaf_sum = 0.0f;

    for (BinaryTiling& tiling : this->tilings)
    {
        tiling.total_power = 0.0f;

        for (int y = 0; y < td.y; ++y)
        {
            float row_mean = 0.0f;
            for (int x = 0; x < td.x; ++x)
            {
                BTTracker tracker;

                auto bounds = tiling.base_tile_bounds(x, y);
                tracker.set_boundaries(bounds.first, bounds.second);

                int i = (y * td.x) + x;
                BinaryTile& tile = tiling.tiles[i];
                float base_power = tiling.recurse_statistics(tile, tracker, this->leaf_sum);

                float a = tile.area(tracker, true);
                row_mean += base_power * a;
                tiling.tile_areas[i] = a;
            }

            tiling.total_power += row_mean;
            tiling.row_means[y] = (row_mean / td.x);
        }

        tiling.total_power /= (td.x * td.y);
    }
}

void BinaryTileCoding::postprocess(bool is_last)
{
    build();

    if (is_last)
    {
        for (BinaryTiling& tiling : this->tilings)
        {
            tiling.tiles.shrink_to_fit();
        }
    }
}

Sample BinaryTileCoding::sample(Pair<float>& pos)
{
    if (this->leaf_sum <= 0.0f)
    {
        Pair<float> random = BinaryTileCoding::random.next2D();

        Sample empty;
        empty.phi = 2.0f * M_PI * random.x;
        empty.theta = std::acos(1.0f - 2.0f * random.y);
        empty.value = 0;
        empty.p = INV_FOURPI;

        return empty;
    }

    int tiling_count = this->tilings.size();

    // [1] Pick one of the tilings with equal weight
    float random = BinaryTileCoding::random.next1D();
    int i = (int)(random * tiling_count);
    BinaryTiling& tiling = this->tilings[i];

    // [2] If there is more than one base tile, use the passed-in (ideally random) 2D position to fetch the right one
    BinaryTile* curr_tile = &tiling.tiles[0];
    Pair<float> x_bounds = tiling.x_bounds;
    Pair<float> y_bounds = tiling.y_bounds;

    if (TILE_DIMS.x * TILE_DIMS.y > 1)
    {
        auto bt_pos = tiling.base_tile_pos(pos);
        curr_tile = &tiling.tiles[(bt_pos.y * TILE_DIMS.x) + bt_pos.x];

        auto bounds = tiling.base_tile_bounds(bt_pos.x, bt_pos.y);
        x_bounds = bounds.first;
        y_bounds = bounds.second;
    }

    // [3] Generate random 1D sample and go deeper as long as the tile isn't a leaf
    BTTracker tracker;
    tracker.set_boundaries(x_bounds, y_bounds);

    while (!curr_tile->is_leaf())
    {
        InternalNode& parent = curr_tile->internal;
        BinaryTile* first = &tiling.tiles[parent.children[0]];
        BinaryTile* second = &tiling.tiles[parent.children[1]];

        if (first->sum == 0 && second->sum == 0)
        {
            break;
        }

        bool h_split = (parent.split_direction == Horizontal);
        Pair<float>& bounds = h_split ? x_bounds : y_bounds;
        float halved = (bounds.x + bounds.y) * 0.5f;

        float power_first = std::max(TC_Epsilon, parent.power[0]);
        float power_second = std::max(TC_Epsilon, parent.power[1]);

        float random = BinaryTileCoding::random.next1D();
        float split = power_first / (power_first + power_second);

        if (random < split)
        {
            bounds.y = halved;
            curr_tile = first;
            tracker.sum = power_first;
        }
        else
        {
            bounds.x = halved;
            curr_tile = second;
            tracker.sum = power_second;
        }

        tracker.increment();
        tracker.set_boundaries(x_bounds, y_bounds);
    }

    // [4] Generate random sample once in a leaf tile
    Pair<float> random2d = BinaryTileCoding::random.next2D();

    x_bounds = tracker.clamped(Horizontal);
    y_bounds = tracker.clamped(Vertical);
    y_bounds.x = std::cos(tiling.transform(y_bounds.x) * M_PI);
    y_bounds.y = std::cos(tiling.transform(y_bounds.y) * M_PI);

    Pair<float> coords(
        x_bounds.x + random2d.x * (x_bounds.y - x_bounds.x),
        y_bounds.x + random2d.y * (y_bounds.y - y_bounds.x)
    );
    coords.y = std::acos(coords.y) * INV_PI;

    coords.x = Util::clamp(coords.x, 0.0f, 1.0f - TC_Epsilon);
    coords.y = Util::clamp(coords.y, 0.0f, 1.0f - TC_Epsilon);

    // [5] Calculate p(x) based on position
    float prob = curr_tile->sum / curr_tile->area(tracker, true);
    for (int ti = 0; ti < tiling_count; ++ti)
    {
        if (ti == i) continue;
        tracker = BTTracker();
        BinaryTile& tile = this->tilings[ti].find_tile(coords, tracker, false);
        prob += tile.sum / tile.area(tracker, true);
    }

    Sample sample;
    sample.phi = coords.x * 2 * M_PI;
    sample.theta = coords.y * M_PI;
    sample.value = 0;
    sample.p = std::max(TC_Epsilon, prob / this->leaf_sum);

    return sample;
}

float BinaryTileCoding::eval(Pair<float>& pos)
{
    if (this->leaf_sum <= 0.0f)
    {
        return INV_FOURPI;
    }

    float prob = 0.0f;
    for (auto& tiling : this->tilings)
    {
        BTTracker tracker;
        BinaryTile& tile = tiling.find_tile(pos, tracker, false);
        prob += tile.sum / tile.area(tracker, true);
    }

    return std::max(TC_Epsilon, (prob / this->leaf_sum));
}

void BinaryTileCoding::wipe()
{
    this->leaf_sum = 0.0f;
    this->split_count = 0;

    for (auto& tiling : this->tilings)
    {
        tiling = BinaryTiling();
    }
}

int BinaryTileCoding::memory()
{
    size_t total_size = sizeof(*this);
    total_size += this->tilings.capacity() * sizeof(BinaryTiling);

    for (const auto& tiling : this->tilings)
    {
        total_size += tiling.tiles.capacity() * sizeof(BinaryTile);
        total_size += tiling.row_means.capacity() * sizeof(float);
        total_size += tiling.tile_areas.capacity() * sizeof(float);
    }

    return total_size;
}

TC_NAMESPACE_END