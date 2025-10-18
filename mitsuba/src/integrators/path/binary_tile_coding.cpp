#include "binary_tile_coding.h"

BTC_NAMESPACE_BEGIN

/* ========== */
/* BinaryTile */
/* ========== */

bool BinaryTile::is_leaf() const
{
    return this->data.tile_type == TileType::Leaf;
}

bool BinaryTile::should_split(const TileTracker& tracker) const
{
    auto sample_count = this->data.sample_count;
    if (sample_count < 2)
    {
        return false;
    }

    // We assume that the population mean is 0.
    float std_err = std::sqrt(var(tracker) / sample_count);

    float t_own = meandev(tracker) / std_err;
    float t_req = TTable::fetch(BinaryTileCoding::CI, sample_count - 1);

    return (t_own > t_req);
}

Direction BinaryTile::split_direction(const TileTracker& tracker) const
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

void BinaryTile::update_statistics(const Sample& sample)
{
    auto samples = this->data.sample_count;
    float value_mean = (samples > 0)
        ? (this->sum / samples)
        : 0;

    auto n = samples + 1;

    float dx = sample.value - value_mean;
    auto dy_x = [&]() { return sample.phi - this->leaf.sample_mean.x; };
    auto dy_y = [&]() { return sample.theta - this->leaf.sample_mean.y; };

    // Update x covariance
    this->leaf.sample_mean.x += dy_x() / n;
    this->leaf.cov.x += dx * dy_x();

    // Update y covariance
    this->leaf.sample_mean.y += dy_y() / n;
    this->leaf.cov.y += dx * dy_y();

    // Update mean deviation
    value_mean += dx / n;
    float dx2 = sample.value - value_mean;
    this->leaf.m2 += dx * dx2;
    this->leaf.diff_sum += std::abs(dx2);
}

void BinaryTile::update_sum(const Sample& sample)
{
    this->sum += sample.value;
    this->data.sample_count++;
}

float BinaryTile::meandev(const TileTracker& tracker) const
{
    if (this->data.sample_count == 0) return 0.0f;

    return area(tracker) * (this->leaf.diff_sum / this->data.sample_count);
}

float BinaryTile::var(const TileTracker& tracker) const
{
    if (this->data.sample_count < 2) return 0.0f;

    float a = area(tracker);
    return a * a * (this->leaf.m2 / (this->data.sample_count - 1));
}

float BinaryTile::mean() const
{
    if (this->data.sample_count == 0)
    {
        return 0;
    }

    return this->sum / this->data.sample_count;
}

float BinaryTile::area(const TileTracker& tracker) const
{
    auto x_clamped = tracker.clamped(Horizontal);
    if (x_clamped.x == x_clamped.y) return 0.0f;

    auto y_clamped = tracker.clamped(Vertical);
    if (y_clamped.x == y_clamped.y) return 0.0f;

    float d_theta = std::cos(M_PI * y_clamped.x) - std::cos(M_PI * y_clamped.y);
    float d_phi = 2 * M_PI * (x_clamped.y - x_clamped.x);

    return d_theta * d_phi;
}

/* ============ */
/* BinaryTiling */
/* ============ */

BinaryTile& BinaryTiling::find_tile(const Pair<float>& pos, bool find_empty)
{
    TileTracker unused;
    return find_tile(pos, unused, find_empty);
}

BinaryTile& BinaryTiling::find_tile(const Pair<float>& pos, TileTracker& tracker, bool find_empty)
{
    const Pair<int> td = BinaryTileCoding::TILE_DIMS;

    // Find right base tile by mapping the pos within the range of the tiling to the range [0, 1]
    Pair<float> pos_warped(pos.x, transform(pos.y, true));
    auto x_warped = Converter::map(pos_warped.x).from({ this->x_bounds.x, this->x_bounds.y }).to({ 0.0, 1.0 });
    auto y_warped = Converter::map(pos_warped.y).from({ this->y_bounds.x, this->y_bounds.y }).to({ 0.0, 1.0 });

    Pair<int> bt_pos(
        x_warped * td.x,
        y_warped * td.y
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
        Direction split_dir = curr_tile->node.split_direction;
        bool h_split = (split_dir == Horizontal);

        Pair<float>& bounds = h_split ? x_bounds : y_bounds;
        float local = h_split ? pos.x : pos.y;
        float split = (bounds.x + bounds.y) * 0.5;

        if (!find_empty && children_empty(*curr_tile))
        {
            break;
        }

        float split_pos = h_split ? split : transform(split);
        bool first_child = local < split_pos;
        BinaryTile* child = &this->tiles[curr_tile->node.children[!first_child]];

        if (first_child)
        {
            bounds.y = split; // we're in the left or upper subtree
        }
        else
        {
            bounds.x = split; // we're in the right or lower subtree
        }

        curr_tile = child;

        tracker.increment();
    }

    tracker.boundaries(x_bounds, y_bounds);
    return *curr_tile;
}

void BinaryTiling::insert(const Sample& sample)
{
    // Find initial tile
    TileTracker tracker;
    Pair<float> uv = Converter::spherical_to_uv(Pair<float>(sample.phi, sample.theta));
    BinaryTile& tile = find_tile(uv, tracker);

    // Store value & update covariance
    tile.update_statistics(sample);
    tile.update_sum(sample);

    // Only split if necessary
    if (tracker.depth > BinaryTileCoding::MAX_DEPTH || !tile.should_split(tracker)) return;

    tile.node.split_direction = tile.split_direction(tracker);
    uint32_t tiles = this->tiles.size();
    tile.node.children = { tiles, tiles + 1 };
    tile.node.power = { 0.0f, 0.0f };
    tile.data.tile_type = TileType::Node;

    BinaryTile child(Leaf);
    this->tiles.push_back(child);
    this->tiles.push_back(child);
}

float BinaryTiling::recurse_statistics(BinaryTile& curr_tile, TileTracker tracker, float& leaf_sum)
{
    if (curr_tile.is_leaf() || children_empty(curr_tile))
    {
        tracker.y_bounds.x = transform(tracker.y_bounds.x);
        tracker.y_bounds.y = transform(tracker.y_bounds.y);

        float power = curr_tile.mean() * curr_tile.area(tracker);
        leaf_sum += power;

        return power;
    }

    for (int i = 0; i < 2; ++i)
    {
        BinaryTile& child = this->tiles[curr_tile.node.children[i]];
        if (child.mean() == 0) continue;

        TileTracker c_tracker = tracker;
        c_tracker.increment();

        Pair<float>& bounds = (curr_tile.node.split_direction == Horizontal)
            ? c_tracker.x_bounds
            : c_tracker.y_bounds;

        float& value = (i == 0) ? bounds.y : bounds.x;
        value = (bounds.x + bounds.y) * 0.5;

        curr_tile.node.power[i] += recurse_statistics(child, c_tracker, leaf_sum);
    }

    return curr_tile.node.power[0] + curr_tile.node.power[1];
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
    float sum_y = 0.0f; size_t y = 0;
    size_t y_len = BinaryTileCoding::TILE_DIMS.y;
    for (; y < y_len; ++y)
    {
        sum_y += this->row_means[y];
        if (sum_y / (y_len * this->total_base_mean) >= pos.y) break;
    }
    if (y == y_len) y -= 1;

    // Conditional sampling (x-dir)
    float sum_x = 0.0f; size_t x = 0;
    size_t x_len = BinaryTileCoding::TILE_DIMS.x;
    for (; x < x_len; ++x)
    {
        int i = (y * x_len) + x;
        sum_x += this->tiles[i].mean() * this->tile_areas[i];
        if (sum_x / (x_len * this->row_means[y]) >= pos.x) break;
    }
    if (x == x_len) x -= 1;

    return Pair<int>(x, y);
}

inline bool BinaryTiling::children_empty(BinaryTile& tile) const
{
    if (tile.is_leaf()) return false;

    return (this->tiles[tile.node.children[0]].mean() == 0 && this->tiles[tile.node.children[1]].mean() == 0);
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
int BinaryTileCoding::MAX_DEPTH = 10;
TTable::CI BinaryTileCoding::CI = TTable::CI::P950;
Transformation::Domain BinaryTileCoding::TRANSFORM_MODE = Transformation::Domain::Spherical;

void BinaryTileCoding::construct(int tilings, Pair<int> tiles, int max_depth, int eagerness, Transformation::Domain trafo)
{
    assert(tilings > 0);
    assert(tiles.x > 0 && tiles.y > 0);
    assert(max_depth > 0);
    assert(-1 < eagerness && eagerness < 5);

    BinaryTileCoding::TILE_DIMS = tiles;
    BinaryTileCoding::MAX_DEPTH = max_depth;
    BinaryTileCoding::CI = static_cast<TTable::CI>(eagerness);
    BinaryTileCoding::TRANSFORM_MODE = trafo;

    this->tilings = std::vector<BinaryTiling>(tilings);
}

void BinaryTileCoding::preprocess()
{
    // Allocate tiles space
    int base_tiles = TILE_DIMS.x * TILE_DIMS.y;
    int max_cap = (1 << (MAX_DEPTH + 1)) - 1;

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
    auto num_tilings = this->tilings.size();

    float overhead = 0.0;
    if (num_tilings > 1)
    {
        overhead = 1.0 / num_tilings;
    }

    // Store start and end thresholds in each tiling
    for (size_t i = 0; i < num_tilings; ++i)
    {
        BinaryTiling& tiling = this->tilings[i];

        Pair<float> x_range(0 - (i * overhead), 1 + ((num_tilings - 1 - i) * overhead));
        Pair<float> y_range(0 - ((num_tilings - 1 - i) * overhead), 1 + (i * overhead));

        tiling.x_bounds = x_range;
        tiling.y_bounds = y_range;
    }
}

void BinaryTileCoding::store(std::vector<Sample>& samples)
{
#if BTC_MULTITHREADING
    auto insert = [&](BinaryTiling& tiling) {
        for (auto& sample : samples)
        {
            tiling.insert(sample);
        }
        };

    std::vector<std::thread> pool;
    for (BinaryTiling& tiling : this->tilings)
    {
        std::thread t(insert, std::ref(tiling));
        pool.push_back(std::move(t));
    }

    for (auto& t : pool)
    {
        if (t.joinable())
        {
            t.join();
        }
    }
#else
    for (auto& sample : samples)
    {
        for (auto& tiling : this->tilings)
        {
            tiling.insert(sample);
        }
    }
#endif
}

void BinaryTileCoding::postprocess()
{
    const Pair<int> td = BinaryTileCoding::TILE_DIMS;

    for (BinaryTiling& tiling : this->tilings)
    {
        for (int y = 0; y < td.y; ++y)
        {
            float row_mean = 0.0f;
            for (int x = 0; x < td.x; ++x)
            {
                TileTracker tracker;

                auto bounds = tiling.base_tile_bounds(x, y);
                tracker.boundaries(bounds.first, bounds.second);

                int i = (y * td.x) + x;
                BinaryTile& tile = tiling.tiles[i];
                tiling.recurse_statistics(tile, tracker, this->leaf_sum);

                tracker.y_bounds.x = tiling.transform(tracker.y_bounds.x);
                tracker.y_bounds.y = tiling.transform(tracker.y_bounds.y);

                float mu = tile.mean();
                float a = tile.area(tracker);
                row_mean += mu * a;

                tiling.tile_areas[i] = a;
            }

            tiling.total_base_mean += row_mean;
            tiling.row_means[y] = (row_mean / td.x);
        }

        tiling.total_base_mean /= (td.x * td.y);
        tiling.tiles.shrink_to_fit();
    }
}

Sample BinaryTileCoding::sample(Pair<float>& pos)
{
    int tiling_count = this->tilings.size();

    // [1] Pick one of the tilings with equal weight
    float random = BinaryTileCoding::random.next1D();
    int i = random * tiling_count;
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
    TileTracker tracker;
    tracker.boundaries(x_bounds, y_bounds);

    while (!curr_tile->is_leaf())
    {
        BinaryTile* first = &tiling.tiles[curr_tile->node.children[0]];
        BinaryTile* second = &tiling.tiles[curr_tile->node.children[1]];

        if (first->mean() == 0 && second->mean() == 0)
        {
            break;
        }

        bool h_split = (curr_tile->node.split_direction == Horizontal);
        Pair<float>& bounds = h_split ? x_bounds : y_bounds;
        float halved = (bounds.x + bounds.y) * 0.5;

        float power_first = std::max(BTC_Epsilon, curr_tile->node.power[0]);
        float power_second = std::max(BTC_Epsilon, curr_tile->node.power[1]);

        float random = BinaryTileCoding::random.next1D();
        float split = power_first / (power_first + power_second);

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

        tracker.increment();
        tracker.boundaries(x_bounds, y_bounds);
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

    coords.x = Util::clamp(coords.x, 0.0f, 1.0f - BTC_Epsilon);
    coords.y = Util::clamp(coords.y, 0.0f, 1.0f - BTC_Epsilon);

    // [5] Calculate p(x) based on position
    float prob = curr_tile->mean();
    for (int ti = 0; ti < tiling_count; ++ti)
    {
        if (ti == i) continue;
        BinaryTile& tile = this->tilings[ti].find_tile(coords, false);
        prob += tile.mean();
    }
    prob = std::max(BTC_Epsilon, (prob / this->leaf_sum));

    Sample sample;
    sample.p = prob;
    sample.phi = coords.y * M_PI;
    sample.theta = coords.x * 2 * M_PI;

    return sample;
}

float BinaryTileCoding::eval(Pair<float>& pos)
{
    float prob = 0.0f;
    for (auto& tiling : this->tilings)
    {
        BinaryTile& tile = tiling.find_tile(pos, false);
        prob += tile.mean();
    }

    return std::max(BTC_Epsilon, (prob / this->leaf_sum));
}

void BinaryTileCoding::wipe()
{
    this->leaf_sum = 0.0f;
    for (auto& tiling : this->tilings)
    {
        tiling = BinaryTiling();
    }
}

int BinaryTileCoding::memory()
{
    size_t size_self = sizeof(*this); // base layer size
    size_t size_tilings = this->tilings.capacity() * sizeof(BinaryTiling); // #tilings * base tiling size

    size_t size_tiles = 0;
    for (const auto& tiling : this->tilings)
    {
        // #tiles * constant tile size
        size_tiles += tiling.tiles.capacity() * sizeof(BinaryTile);
    }

    return size_self + size_tilings + size_tiles;
}

BTC_NAMESPACE_END