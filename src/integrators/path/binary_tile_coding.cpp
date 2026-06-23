#include "binary_tile_coding.h"

TC_NAMESPACE_BEGIN

/* ========== */
/* BinaryTile */
/* ========== */

void BinaryTile::update(const Sample& sample)
{
    if (UNLIKELY(!is_leaf())) return;

    auto samples = this->data.sample_count;
    float value_mean = (samples > 0) ? (this->sum / samples) : 0;
    const float inv_pdf = 1.0f / sample.pdf;
    const float contribution = sample.value * inv_pdf;
    auto n = samples + 1;
    const float inv_n = 1.0f / (float)n;

    float dx = contribution - value_mean;
    auto dy_x = [&]() { return sample.phi - this->leaf.sample_mean.x; };
    auto dy_y = [&]() { return sample.theta - this->leaf.sample_mean.y; };

    // Update covariance
    this->leaf.sample_mean.x += dy_x() * inv_n;
    this->leaf.cov.x += dx * dy_x(); // * std::sin(sample.theta);

    this->leaf.sample_mean.y += dy_y() * inv_n;
    this->leaf.cov.y += dx * dy_y();

    // Update mean deviation
    value_mean += dx * inv_n;
    float dx2 = contribution - value_mean;
    float diff_squared = dx * dx2;
    this->leaf.m2 += diff_squared;
    this->leaf.ad_sum += std::sqrt(std::max(0.0f, diff_squared));

    this->data.sample_count++;
    this->sum += contribution;
}

Direction BinaryTile::split_direction(const BTTracker& tracker) const
{
    auto x_norm = tracker.clamped(Horizontal);
    auto y_norm = tracker.clamped(Vertical);

    const float covar_x = (x_norm.x != x_norm.y) ? std::abs(this->leaf.cov.x) : 0.0f;
    const float covar_y = (y_norm.x != y_norm.y) ? std::abs(this->leaf.cov.y) : 0.0f;

    return (covar_x > covar_y) ? Horizontal : Vertical;
}

bool BinaryTile::should_split() const
{
    auto sample_count = this->data.sample_count;
    if (sample_count < 2) return false;

    float mad = meandev();
    if (mad <= 0.0f) return false;

    float ss_y = this->leaf.m2 - (sample_count * mad * mad);
    ss_y = std::max(0.0f, ss_y);

    float var_mad = ss_y / (sample_count - 1);
    if (var_mad <= 0.0f) return true;

    float std_err_squared = var_mad / sample_count;
    float t_req = Lookup::TTable::fetch(BinaryTileCoding::CI, sample_count - 1);

    return ((mad * mad) / std_err_squared > (t_req * t_req));
}

float BinaryTile::meandev() const {
    return (this->data.sample_count == 0) ? 0.0f : (this->leaf.ad_sum / this->data.sample_count);
}

float BinaryTile::var() const {
    return (this->data.sample_count < 2) ? 0.0f : (this->leaf.m2 / (this->data.sample_count - 1));
}

float BinaryTile::area(const BTTracker& tracker, bool transform)
{
    auto x_clamped = tracker.clamped(Horizontal);
    float dx = x_clamped.y - x_clamped.x;
    if (dx == 0.0f) return 0.0f;

    auto y_clamped = tracker.clamped(Vertical);
    float dy = y_clamped.y - y_clamped.x;
    if (dy == 0.0f) return 0.0f;

    double d_theta;

    if (!transform)
    {
        d_theta = std::cos(y_clamped.x * M_PI) - std::cos(y_clamped.y * M_PI);
        double d_phi = 2.0 * M_PI * dx;

        return (float)(d_theta * d_phi);
    }

    if (LIKELY(BinaryTileCoding::TRANSFORM_MODE == DomainTransform::Spherical))
    {
        d_theta = 2.0 * dy; // Math cancels out, resulting in a purely linear delta.
    }
    else if (BinaryTileCoding::TRANSFORM_MODE == DomainTransform::Cosine) // Cosine mapping is weighted, so we still must evaluate the trig
    {
        double theta_min = DomainTransform::apply<float>(DomainTransform::Cosine, y_clamped.x, false) * M_PI;
        double theta_max = DomainTransform::apply<float>(DomainTransform::Cosine, y_clamped.y, false) * M_PI;
        d_theta = std::cos(theta_min) - std::cos(theta_max);
    }
    else
    {
        d_theta = dy; // Planar fallback
    }

    double d_phi = 2.0 * M_PI * dx;
    return (float)(d_theta * d_phi);
}

/* ============ */
/* BinaryTiling */
/* ============ */

BinaryTile& BinaryTiling::find_tile(const Point2& pos, bool find_empty) {
    BTTracker unused;
    return find_tile(pos, unused, find_empty);
}

BinaryTile& BinaryTiling::find_tile(const Point2& pos, BTTracker& tracker, bool find_empty)
{
    const Point2i td = BinaryTileCoding::TILE_DIMS;

    // Map Geometric -> Uniform domain
    Point2 pos_warped(pos.x, DomainTransform::apply<float>(BinaryTileCoding::TRANSFORM_MODE, pos.y, true));
    auto x_warped = Util::map(pos_warped.x).from({ this->x_bounds.x, this->x_bounds.y }).to({ 0.0f, 1.0f });
    auto y_warped = Util::map(pos_warped.y).from({ this->y_bounds.x, this->y_bounds.y }).to({ 0.0f, 1.0f });

    Point2i bt_pos(x_warped * td.x, y_warped * td.y);
    int i = (bt_pos.y * td.x) + bt_pos.x;
    BinaryTile* curr_tile = &this->tiles[i];

    auto bounds = base_tile_bounds(bt_pos.x, bt_pos.y);
    Point2 x_bounds = bounds.first;
    Point2 y_bounds = bounds.second;

    while (!curr_tile->is_leaf())
    {
        InternalNode& parent = curr_tile->internal;
        Direction split_dir = (Direction) curr_tile->data.split_direction;
        bool h_split = (split_dir == Horizontal);

        Point2& c_bounds = h_split ? x_bounds : y_bounds;
        Float local = h_split ? pos.x : pos.y;
        Float split = (c_bounds.x + c_bounds.y) * 0.5f;

        if (!find_empty && children_empty(*curr_tile)) break;

        Float split_pos = h_split ? split : DomainTransform::apply<float>(BinaryTileCoding::TRANSFORM_MODE, split);
        uint32_t child_idx = (local >= split_pos);
        reinterpret_cast<Float*>(&c_bounds)[1 - child_idx] = split;

        tracker.increment();
        tracker.sum = parent.power[child_idx];
        curr_tile = &this->tiles[parent.children[child_idx]];
    }

    tracker.set_boundaries(x_bounds, y_bounds);

    if (tracker.depth == 0)
    {
        tracker.sum = curr_tile->sum;
    }

    return *curr_tile;
}

void BinaryTiling::insert(const Sample& sample, std::atomic<uint32_t>& split_count, float leaf_sum)
{
    // Find initial tile
    BTTracker tracker;
    Point2 uv = Util::spherical_to_uv(Point2(sample.phi, sample.theta));
    BinaryTile& tile = find_tile(uv, tracker);

    // Store value & update covariance
    tile.update(sample);

    // Only split if possible & necessary
    bool below_max_splits = (split_count.load(std::memory_order_relaxed) < BinaryTileCoding::MAX_SPLITS);
    bool split_needed = tile.should_split();

    if (!(below_max_splits && split_needed)) return;

    Direction split_dir = tile.split_direction(tracker);
    float width = 0.0f;
    float cov = 0.0f;

    if (split_dir == Horizontal)
    {
        float delta = tracker.clamped(Horizontal).y - tracker.clamped(Horizontal).x;
        width = delta * 2.0f * M_PI;
        cov = tile.leaf.cov.x;
    }
    else
    {
        float delta = tracker.clamped(Vertical).y - tracker.clamped(Vertical).x;
        width = delta * M_PI;
        cov = tile.leaf.cov.y;
    }

    // This section here is still highly experimental. Please just comment the entire if block (but leave the R at 0.5) if your output causes issues.
    float R = 0.5f;
    if (tile.sum > 0.0f && width > 0.0f && tile.data.sample_count > 1) 
    {
        float N = (float) tile.data.sample_count;
        float MD = tile.meandev();
        float cov_norm = std::abs(cov) / width; 
        
        float var_expl = 4.0f * MD * (cov_norm / N);
        float var_total = tile.leaf.m2 / (N - 1.0f);
        
        float rho = 0.0f;
        if (var_total > 0.0f)
        {
            rho = std::max(0.0f, std::min(var_expl / var_total, 1.0f));
        }

        float sign_cov = (cov > 0.0f) ? 1.0f : ((cov < 0.0f) ? -1.0f : 0.0f);
        float R_raw = 0.5f - sign_cov * (std::sqrt(N * MD * cov_norm) / tile.sum);
        
        R = 0.5f + rho * (R_raw - 0.5f);
        R = std::max(0.05f, std::min(R, 0.95f)); // Clamp for safety
    }

    LeafNode old_leaf = tile.leaf;
    float old_sum = tile.sum;
    uint32_t old_count = tile.data.sample_count;

    // Create new children based on ratio/slope
    BinaryTile child0(Leaf);
    child0.sum = old_sum * R;
    child0.data.sample_count = (uint32_t) old_count * R;
    child0.leaf.m2 = old_leaf.m2 * R;
    child0.leaf.ad_sum = old_leaf.ad_sum * R;
    child0.leaf.sample_mean = old_leaf.sample_mean;

    BinaryTile child1(Leaf);
    child1.sum = old_sum * (1.0f - R);
    child1.data.sample_count = old_count - child0.data.sample_count;
    child1.leaf.m2 = old_leaf.m2 * (1.0f - R);
    child1.leaf.ad_sum = old_leaf.ad_sum * (1.0f - R);
    child1.leaf.sample_mean = old_leaf.sample_mean;

    float mean_shift = width * 0.25f;
    if (split_dir == Horizontal)
    {
        child0.leaf.sample_mean.x -= mean_shift;
        child1.leaf.sample_mean.x += mean_shift;
    } 
    else
    {
        child0.leaf.sample_mean.y -= mean_shift;
        child1.leaf.sample_mean.y += mean_shift;
    }

    // adds memory for the children on demand instead of using memory from the 
    // worst-case preallocation.
    // now only appends the two children when they are needed. 
    // this results in a dependency where the parent tile must be written beofre push_back
    // because push_back can re-allocate tiles and invalidate the reference 

    /* origional code 
    
    uint32_t child_idx = this->next_node_idx.fetch_add(2, std::memory_order_relaxed);
    if (child_idx + 1 < this->tiles.size())
    {
        this->tiles[child_idx] = child0;
        this->tiles[child_idx + 1] = child1;
        InternalNode parent;
        parent.children = { child_idx, child_idx + 1 };
        tile.internal = parent;
        tile.data.split_direction = split_dir;
        tile.data.tile_type = TileType::Internal;
        split_count.fetch_add(1, std::memory_order_relaxed);
    }
    */

    uint32_t child_idx = (uint32_t) this->tiles.size();

    InternalNode parent;
    parent.children = { child_idx, child_idx + 1 };

    tile.internal = parent;
    tile.data.split_direction = split_dir;
    tile.data.tile_type = TileType::Internal;

    this->tiles.push_back(child0);
    this->tiles.push_back(child1);
    this->next_node_idx.store(child_idx + 2, std::memory_order_relaxed);

    split_count.fetch_add(1, std::memory_order_relaxed);
}

float BinaryTiling::recurse_statistics(BinaryTile& curr_tile, BTTracker tracker, float& leaf_sum)
{
    if (curr_tile.is_leaf() || children_empty(curr_tile))
    {
        float power = curr_tile.sum;
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

        Point2& bounds = (curr_tile.data.split_direction == Horizontal)
            ? c_tracker.x_bounds
            : c_tracker.y_bounds;

        Float& value = reinterpret_cast<Float*>(&bounds)[1 - i];
        value = (bounds.x + bounds.y) * 0.5f;
        
        parent.power[i] += recurse_statistics(child, c_tracker, leaf_sum);
    }

    curr_tile.sum = parent.power[0] + parent.power[1];
    return curr_tile.sum;
}

std::pair<Point2, Point2> BinaryTiling::base_tile_bounds(int x, int y) const
{
    Point2i td = BinaryTileCoding::TILE_DIMS;
    float x_step = (this->x_bounds.y - this->x_bounds.x) / td.x;
    float y_step = (this->y_bounds.y - this->y_bounds.x) / td.y;

    Point2 x_bounds(this->x_bounds.x + (x * x_step), this->x_bounds.x + ((x + 1) * x_step));
    Point2 y_bounds(this->y_bounds.x + (y * y_step), this->y_bounds.x + ((y + 1) * y_step));

    return { x_bounds, y_bounds };
}

Point2i BinaryTiling::base_tile_pos(const Point2& pos)
{
    // Marginal sampling (y-dir)
    size_t y = 0;
    float sum_y = 0.0f;
    size_t y_len = BinaryTileCoding::TILE_DIMS.y;

    for (; y < y_len; ++y) {
        sum_y += this->row_sums[y];
        if (sum_y / this->total_power >= pos.y) break;
    }
    if (y == y_len) y -= 1;

    // Conditional sampling (x-dir)
    size_t x = 0;
    float sum_x = 0.0f;
    size_t x_len = BinaryTileCoding::TILE_DIMS.x;

    for (; x < x_len; ++x) {
        int i = (y * x_len) + x;
        float power = this->tiles[i].sum;
        sum_x += power;

        if (sum_x / this->row_sums[y] >= pos.x) break;
    }
    if (x == x_len) x -= 1;

    return Point2i(x, y);
}

inline bool BinaryTiling::children_empty(BinaryTile& tile) const {
    if (tile.is_leaf()) return false;
    return (this->tiles[tile.internal.children[0]].sum == 0) && (this->tiles[tile.internal.children[1]].sum == 0);
}

/* ================ */
/* BinaryTileCoding */
/* ================ */

Point2i BinaryTileCoding::TILE_DIMS                         = Point2i(1, 1);
float BinaryTileCoding::TILING_EXCESS                       = 0.1f;
uint32_t BinaryTileCoding::MAX_SPLITS                       = UINT32_MAX;
Lookup::CI BinaryTileCoding::CI                             = Lookup::CI::P999;
DomainTransform::Type BinaryTileCoding::TRANSFORM_MODE      = DomainTransform::Spherical;

void BinaryTileCoding::construct(const BTCArguments& args)
{
    TC_Assert(args.tilings > 0);
    TC_Assert(args.tiles.x > 0 && args.tiles.y > 0);
    TC_Assert(args.eagerness > -1 && args.eagerness < 5);

    BinaryTileCoding::TILE_DIMS = args.tiles;
    BinaryTileCoding::TILING_EXCESS = args.excess;
    BinaryTileCoding::MAX_SPLITS = args.max_splits;
    BinaryTileCoding::CI = static_cast<Lookup::CI>(args.eagerness);
    BinaryTileCoding::TRANSFORM_MODE = args.transformation;

    this->tilings = std::vector<BinaryTiling>(args.tilings);
}

void BinaryTileCoding::preprocess()
{
    int base_tiles = TILE_DIMS.x * TILE_DIMS.y;

        // origional code
        // std::size_t total_expected_tiles = base_tiles + (MAX_SPLITS * 2); 
        
    for (auto& tiling : this->tilings)
    {
        tiling.row_sums.resize(TILE_DIMS.y, 0.0f);

        // origional code
        // tiling.tiles.resize(total_expected_tiles, BinaryTile(Leaf));
        
        // The child tiles are now added when needed in the ::insert() function
        // this used to preallocate the worst case memory usage per tiling
        // so spatial leafs cost ~1000kb per depending on settings, resulting in exessive memory usage
        tiling.tiles.resize(base_tiles, BinaryTile(Leaf));
        tiling.next_node_idx.store(base_tiles, std::memory_order_relaxed);
    }

    auto tiling_count = this->tilings.size();
    float shift = (tiling_count > 1) ? (TILING_EXCESS / (tiling_count - 1)) : 0.0f;

    for (size_t i = 0; i < tiling_count; ++i)
    {
        BinaryTiling& tiling = this->tilings[i];
        tiling.x_bounds = Point2(0 - (i * shift), 1 + ((tiling_count - 1 - i) * shift)); 
        tiling.y_bounds = Point2(0 - ((tiling_count - 1 - i) * shift), 1 + (i * shift));
    }
}

void BinaryTileCoding::store(Sample& sample)
{
    for (auto& tiling : this->tilings)
    {
        tiling.insert(sample, this->split_count, this->leaf_sum);
    }
}

void BinaryTileCoding::build()
{
    const Point2i td = BinaryTileCoding::TILE_DIMS;
    this->leaf_sum = 0.0f;

    #pragma omp parallel for schedule(dynamic)
    for (size_t t = 0; t < this->tilings.size(); ++t)
    {
        BinaryTiling& tiling = this->tilings[t];
        tiling.total_power = 0.0f;
        float local_leaf_sum = 0.0f;

        for (int y = 0; y < td.y; ++y)
        {
            float row_sum = 0.0f;
            for (int x = 0; x < td.x; ++x)
            {
                BTTracker tracker;

                auto bounds = tiling.base_tile_bounds(x, y);
                tracker.set_boundaries(bounds.first, bounds.second);

                int i = (y * td.x) + x;
                BinaryTile& tile = tiling.tiles[i];
                float base_power = tiling.recurse_statistics(tile, tracker, local_leaf_sum);

                row_sum += base_power;
            }

            tiling.total_power += row_sum;
            tiling.row_sums[y] = row_sum;
        }

        #pragma omp atomic
        this->leaf_sum += local_leaf_sum;
    }
}

void BinaryTileCoding::postprocess(bool last_iteration)
{
    build();

    if (last_iteration)
    {
        #pragma omp parallel for schedule(dynamic)
        for (size_t t = 0; t < this->tilings.size(); ++t)
        {
            this->tilings[t].tiles.shrink_to_fit();
        }
    }
}

Sample BinaryTileCoding::sample(Point2 pos)
{
    if (this->leaf_sum <= 0.0f)
    {
        Point2 random = get_thread_random().next2D();
        Sample empty = {
            .phi = 2.0f * (float) M_PI * random.x,
            .theta = std::acos(1.0f - 2.0f * random.y),
            .value = 0.0f,
            .pdf = INV_FOURPI * std::sin(random.y)
        };
        return empty;
    }

    int tiling_count = this->tilings.size();

    // [1] Pick one of the tilings with equal weight
    Float random = get_thread_random().next1D();
    int i = random * tiling_count;
    BinaryTiling& tiling = this->tilings[i];

    // [2] Fetch base tile
    BinaryTile* curr_tile = &tiling.tiles[0];
    Point2 x_bounds = tiling.x_bounds;
    Point2 y_bounds = tiling.y_bounds;
    
    if (TILE_DIMS.x * TILE_DIMS.y > 1)
    {
        auto bt_pos = tiling.base_tile_pos(pos);
        curr_tile = &tiling.tiles[(bt_pos.y * TILE_DIMS.x) + bt_pos.x];

        auto bounds = tiling.base_tile_bounds(bt_pos.x, bt_pos.y);
        x_bounds = bounds.first;
        y_bounds = bounds.second;
    }

    // [3] Tree traversal
    BTTracker tracker;
    tracker.set_boundaries(x_bounds, y_bounds);
    tracker.sum = curr_tile->sum;

    while (!curr_tile->is_leaf())
    {
        InternalNode& parent = curr_tile->internal;
        float power_first = parent.power[0];
        float power_second = parent.power[1];

        if (power_first == 0 && power_second == 0)
        {
            break;
        }

        bool h_split = (curr_tile->data.split_direction == Horizontal);
        Point2& bounds = h_split ? x_bounds : y_bounds;
        Float halved = (bounds.x + bounds.y) * 0.5;

        Float random = get_thread_random().next1D();
        float total_power = power_first + power_second;
        uint32_t child_idx = ((random * total_power) >= power_first);

        reinterpret_cast<Float*>(&bounds)[1 - child_idx] = halved;
        tracker.sum = parent.power[child_idx];
        curr_tile = &tiling.tiles[parent.children[child_idx]];

        tracker.increment();
        tracker.set_boundaries(x_bounds, y_bounds);
    }

    // [4] Generate random sample in leaf tile
    Point2 random2d = get_thread_random().next2D();

    x_bounds = tracker.clamped(Horizontal);
    y_bounds = tracker.clamped(Vertical);
    y_bounds.x = std::cos(DomainTransform::apply<float>(TRANSFORM_MODE, y_bounds.x) * M_PI);
    y_bounds.y = std::cos(DomainTransform::apply<float>(TRANSFORM_MODE, y_bounds.y) * M_PI);

    Point2 coords(
        x_bounds.x + random2d.x * (x_bounds.y - x_bounds.x),
        y_bounds.x + random2d.y * (y_bounds.y - y_bounds.x)
    );

    coords.y = Util::clamp(coords.y, -1.0f, 1.0f);
    coords.y = std::acos(coords.y) * INV_PI;

    // [5] Calculate PDF based on area
    float prob = tracker.sum / curr_tile->area(tracker, true);
    for (int ti = 0; ti < tiling_count; ++ti)
    {
        if (ti == i) continue;
        tracker = BTTracker();
        BinaryTile& tile = this->tilings[ti].find_tile(coords, tracker, false);
        prob += tracker.sum / tile.area(tracker, true);
    }
    prob = std::max(TC_Epsilon, prob / this->leaf_sum);

    Sample sample = {
        .phi = (coords.x * 2.0f * (float) M_PI),
        .theta = (coords.y * (float) M_PI),
        .value = 0.0f,
        .pdf = prob
    };
    return sample;
}

Float BinaryTileCoding::eval(Point2 pos)
{
    if (this->leaf_sum <= 0.0f) return INV_FOURPI;

    float prob = 0.0f;
    for (auto& tiling : this->tilings)
    {
        BTTracker tracker;
        BinaryTile& tile = tiling.find_tile(pos, tracker, false);
        prob += tracker.sum / tile.area(tracker, true);
    }

    return std::max(TC_Epsilon, (prob / this->leaf_sum));
}

void BinaryTileCoding::wipe()
{
    this->leaf_sum = 0.0f;
    this->split_count = 0;

    int base_tiles = TILE_DIMS.x * TILE_DIMS.y;

    for (auto& tiling : this->tilings)
    {
        for (int i = 0; i < base_tiles; ++i) {
            tiling.tiles[i] = BinaryTile(Leaf);
        }

        std::fill(tiling.row_sums.begin(), tiling.row_sums.end(), 0.0f);
        
        tiling.total_power = 0.0f;
        tiling.next_node_idx.store(base_tiles, std::memory_order_relaxed);
    }
}

int BinaryTileCoding::memory()
{
    size_t total_size = sizeof(*this);
    total_size += this->tilings.capacity() * sizeof(BinaryTiling);

    for (const auto& tiling : this->tilings)
    {
        total_size += tiling.tiles.capacity() * sizeof(BinaryTile);
        total_size += tiling.row_sums.capacity() * sizeof(float);
    }

    return total_size;
}

TC_NAMESPACE_END