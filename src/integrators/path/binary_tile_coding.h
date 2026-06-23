#pragma once

#include "tc_util.h"

TC_NAMESPACE_BEGIN

// Riley Swinson: 
// I editied portions of this header so I can use atomics, and then I re-define the 
// copy, move, and copy assign constructors so they are compatible with automics

/**
 * @brief Used to specify the direction of an entity, e.g. a split, amongst other things.
 * 
 * In the case of a split:
 * HORIZONTAL = divides horizontally, i.e. the split goes from top to bottom.
 * VERTICAL = divides vertically, i.e. the split goes from left to right.
 */
enum Direction {
    Horizontal,
    Vertical
};

/**
 * @brief Helper struct to track tile data.
 */
struct BTTracker {
    int depth = 0;
    float sum = 0;
    Point2 x_bounds;
    Point2 y_bounds;

    inline void increment()
    {
        this->depth++;
    }

    inline void set_boundaries(const Point2 x, const Point2 y)
    {
        this->x_bounds = x;
        this->y_bounds = y;
    }

    inline Point2 clamped(const Direction dir) const
    {
        Point2 bounds = (dir == Horizontal) ? this->x_bounds : this->y_bounds;
        return Point2(
            Util::clamp(bounds.x, 0.0f, 1.0f),
            Util::clamp(bounds.y, 0.0f, 1.0f)
        );
    }
};

enum TileType {
    Internal,
    Leaf
};

struct BTBitfield {
    uint32_t sample_count : 30;     // Number of samples (max.: 2^30 ~ 1.073bn)
    uint32_t tile_type : 1;         // 0 = Leaf, 1 = Internal
    uint32_t split_direction : 1;   // 0 = Horizontal, 1 = Vertical

    BTBitfield() : sample_count(0), tile_type(Leaf), split_direction(Horizontal) {};
};

/**
 * @brief Tracks the statistics of a tile.
 *
 * Upon sample storage, the statistics are updated via a 1-pass Welfold algorithm and used to determine...
 * (a) if a tile should be split, by conducting a Student's 1-sample t-test.
 * (b) how a tile should be split, by looking at the covariances in both directions (phi, theta).
 */
struct LeafNode {
    // The covariance in x and y direction to determine how to split.
    Point2 cov;
    // Sample mean required for covariance updates.
    Point2 sample_mean;
    // Sum of squared distances from the mean.
    float m2;
    // Sum of absolute deviations used for MD metric.
    float ad_sum;
};

/**
 * @brief Stores relevant data about child tiles.
 */
struct InternalNode {
    // The total radiance in each child (including sub-children).
    std::array<float, 2> power;
    // The indices of the children. Use in combination with the tiles vector.
    std::array<uint32_t, 2> children;
};

/**
 * @brief Tile in a tiling.
 * 
 * The most low-level entity in the 'Binary Tile Coding' data structure, storing radiance information
 * in a specific area of the sample space (and beyond). Can be both leaf and non-leaf ("internal").
 * Normally, a binary tile is initialized as a leaf and collects statistics about samples falling into
 * its region until the should_split() method determines that the data implicates that the tile should
 * be split, after which it turns into an internal node.
 * 
 * To access a child, use the tiles vector in combination with the child indices.
 */
struct alignas(32) BinaryTile {
    // Total radiance that arrived in this tile.
    float sum = 0.0f;
    // Sample count + tile flag
    BTBitfield data = BTBitfield();

    union {
        InternalNode internal;
        LeafNode leaf;
    };

    BinaryTile()
    {
        this->data = BTBitfield();
        this->data.tile_type = TileType::Leaf;
        this->sum = 0.0f;
        
        this->leaf.cov = Point2(0.0f);
        this->leaf.sample_mean = Point2(0.0f);
        this->leaf.m2 = 0.0f;
        this->leaf.ad_sum = 0.0f;
    }

    BinaryTile(const TileType type)
    {
        this->data = BTBitfield();
        this->data.tile_type = type;
        this->sum = 0.0f;

        if (type == TileType::Internal)
        {
            this->internal.power = { 0.0f, 0.0f };
            this->internal.children = { UINT32_MAX, UINT32_MAX };
        }
        else
        {
            this->leaf.cov = Point2(0.0f);
            this->leaf.sample_mean = Point2(0.0f);
            this->leaf.m2 = 0.0f;
            this->leaf.ad_sum = 0.0f;
        }
    }

    /// Correctly updates all required statistics & total radiance.
    inline void update(const Sample& sample);

    /// Returns the split direction of a leaf tile by calculating the absolute covariance in both x (phi) and y (theta) direction.
    Direction split_direction(const BTTracker& tracker) const;

    /// Determines if a leaf tile should be split by performing a one-sample t-test, trying to see if we can reject our null hypothesis that the MD is significantly different from 0.
    bool should_split() const;

    /// Returns the mean deviation of the current leaf.
    float meandev() const;

    /// Returns the sample variance of the current leaf.
    float var() const;

    /// Returns the area of this tile.
    float area(const BTTracker& tracker, bool transform = false);

    /// Returns whether the current tile is a leaf.
    inline bool is_leaf() const { return this->data.tile_type == TileType::Leaf; }
};

/**
 * @brief Tiling containing tiles.
 * 
 * Each 'Binary Tiling' can be understood as a collection of tiles that act as binary trees.
 * The number of binary trees depends on the initialization, and is guaranteed to be x * y, where x is
 * the number of tiles in x-direction and y the #tiles in y-direction. If only one tile is stored during
 * the entire construction process, the tiling will only contain a single tree. If, e.g., a tiling is
 * however initialized with x = y = 4, it will contain 16 b-trees, as each base tile acts as its own tree.
 */
struct BinaryTiling {
    Point2 x_bounds;
    Point2 y_bounds;
    std::vector<BinaryTile> tiles;
    std::atomic<uint32_t> next_node_idx{0};

    /* Information used for marginal and conditional sampling (CDF-based sampling) */
    float total_power = 0.0f;
    std::vector<float> row_sums;

    BinaryTiling() = default;

    // Atomic next_node_idx deletes the implicit copy/move; restore them manually
    // so BinaryTileCoding can deep-copy/move its tilings vector.

    BinaryTiling(BinaryTiling&& other) noexcept :
        x_bounds(other.x_bounds), y_bounds(other.y_bounds),
        tiles(std::move(other.tiles)), total_power(other.total_power),
        row_sums(std::move(other.row_sums))
    {
        next_node_idx.store(other.next_node_idx.load(std::memory_order_relaxed), std::memory_order_relaxed);
    }

    // Needed for vector<BinaryTiling> copy-construction inside BinaryTileCoding's copy constructor.
    BinaryTiling(const BinaryTiling& other) :
        x_bounds(other.x_bounds), y_bounds(other.y_bounds),
        tiles(other.tiles), total_power(other.total_power),
        row_sums(other.row_sums)
    {
        next_node_idx.store(other.next_node_idx.load(std::memory_order_relaxed), std::memory_order_relaxed);
    }

    // Needed for vector<BinaryTiling> copy-assignment in m_sampling = m_building and parent->child BTC copy.
    BinaryTiling& operator=(const BinaryTiling& other) {
        if (this != &other) {
            x_bounds = other.x_bounds;
            y_bounds = other.y_bounds;
            tiles = other.tiles;
            total_power = other.total_power;
            row_sums = other.row_sums;
            next_node_idx.store(other.next_node_idx.load(std::memory_order_relaxed), std::memory_order_relaxed);
        }
        return *this;
    }

    // Symmetry with the move constructor; without it std::move falls back to copy.
    BinaryTiling& operator=(BinaryTiling&& other) noexcept {
        if (this != &other) {
            x_bounds = other.x_bounds;
            y_bounds = other.y_bounds;
            tiles = std::move(other.tiles);
            total_power = other.total_power;
            row_sums = std::move(other.row_sums);
            next_node_idx.store(other.next_node_idx.load(std::memory_order_relaxed), std::memory_order_relaxed);
        }
        return *this;
    }

    /// Finds a tile at a given position within the bounds of this tiling. If find_empty is false, the parent will be returned in case both children are empty.
    BinaryTile& find_tile(const Point2& pos, bool find_empty = true);

    /// Finds a tile at a given position within the bounds of this tiling, but with the option to pass a tracker to obtain further data about the tile. If find_empty is false, the parent will be returned in case both children are empty.
    BinaryTile& find_tile(const Point2& pos, BTTracker& counter, bool find_empty = true);

    /// Stores a sample in a binary tiling.
    void insert(const Sample& sample, std::atomic<uint32_t>& split_count, float leaf_sum);

    /// Utility function to recursively add the accumulated sum from all children to their parent tiles.
    float recurse_statistics(BinaryTile& curr_tile, BTTracker tracker, float& leaf_sum);

    /// Utility function to calculate the planar (!) boundaries of a base tile.
    inline std::pair<Point2, Point2> base_tile_bounds(int x, int y) const;

    /// Utility function to obtain the 2D index of a base tile.
    inline Point2i base_tile_pos(const Point2& pos);

    /// Checks if the children of the passed in tile are both empty, i.e., their weighted sum equals 0.
    inline bool children_empty(BinaryTile& tile) const;
};

/**
 * @brief Uppermost layer of the Binary Tile Coding (BTC).
 */
struct BinaryTileCoding {
    static Lookup::CI                   CI;
    static DomainTransform::Type        TRANSFORM_MODE;
    static uint32_t                     MAX_SPLITS;
    static Point2i                      TILE_DIMS;
    static float                        TILING_EXCESS;

    float leaf_sum = 0.0f;
    std::atomic<uint32_t> split_count{0};

    // Adding any user-declared constructor below would suppress the implicit default constructor;
    // BTCAdapter relies on default-constructing m_building/m_sampling so I'm just keeping it .
    BinaryTileCoding() = default;
    ~BinaryTileCoding() {}

    // Atomic split_count deletes the implicit copy/move; restore them manually so
    // BTCAdapter can do m_sampling = m_building (build) and nodes[i].btc = cur.btc
    // which lets spatial-tree subdivision children inherit the spatial parent's distribution.

    BinaryTileCoding(const BinaryTileCoding& other) :
        leaf_sum(other.leaf_sum), tilings(other.tilings)
    {
        split_count.store(other.split_count.load(std::memory_order_relaxed), std::memory_order_relaxed);
    }

    BinaryTileCoding& operator=(const BinaryTileCoding& other) {
        if (this != &other) {
            leaf_sum = other.leaf_sum;
            tilings = other.tilings;
            split_count.store(other.split_count.load(std::memory_order_relaxed), std::memory_order_relaxed);
        }
        return *this;
    }

    // Move ops so std::move actually steals the tilings vector instead of deep-copying.
    BinaryTileCoding(BinaryTileCoding&& other) noexcept :
        leaf_sum(other.leaf_sum), tilings(std::move(other.tilings))
    {
        split_count.store(other.split_count.load(std::memory_order_relaxed), std::memory_order_relaxed);
    }

    BinaryTileCoding& operator=(BinaryTileCoding&& other) noexcept {
        if (this != &other) {
            leaf_sum = other.leaf_sum;
            tilings = std::move(other.tilings);
            split_count.store(other.split_count.load(std::memory_order_relaxed), std::memory_order_relaxed);
        }
        return *this;
    }

    void construct(const BTCArguments& args);
    void preprocess();
    void store(Sample& sample);
    void postprocess(bool is_last = false);

    Sample sample(Point2 pos);
    float eval(Point2 pos);
    void wipe();
    int memory();

    // Read-only access to the private tilings vector; needed by BTCAdapter::dump()
    // for serialization.
    const std::vector<BinaryTiling>& getTilings() const { return tilings; }

    // Total tile count across all tilings just used for terminal logs 
    size_t nodeCount() const {
        size_t count = 0;
        for (const auto& t : tilings) count += t.tiles.size();
        return count;
    }

private:
    std::vector<BinaryTiling> tilings;
    void build();
};

TC_NAMESPACE_END