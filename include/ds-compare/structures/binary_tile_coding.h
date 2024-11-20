#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_)
#define __DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_

#include <ds-compare/ds.h>

MTS_NAMESPACE_BEGIN

/**
 * @brief Used to specify the split direction of a tile.
 * 
 * HORIZONTAL = horizontal split, split from top to bottom.
 * VERTICAL = vertical split, split from left to right.
 */
enum SplitDirection {
    HORIZONTAL,
    VERTICAL
};

/**
 * @brief Helper struct to track the depth of a tile.
 */
struct DepthCounter {
    int horizontal = 0;
    int vertical = 0;

    inline int depth() const { return this->horizontal + this->vertical; }
    inline void increment(const SplitDirection split) { if (split == HORIZONTAL) this->horizontal++; else this->vertical++; }
};

/**
 * @brief Tile in a tiling.
 * 
 * The most low-level entity in the 'Binary Tile Coding' data structure, storing luminance information
 * in a specific area of the sample space (and beyond). Can be both leaf and non-leaf based on the two
 * sub-leaf indices it may hold. To access a child, use the tiles vector in combination with the child
 * vertices.
 */
struct BinaryTile {
    /* ==== Statistics (24 bytes) ==== */

    Point2f cov = Point2f(0.0f);
    Point2f sample_mean = Point2f(0.0f);
    float m2 = 0;
    float diff_sum = 0;

    /* ==== Data (16 bytes) ==== */

    uint32_t idx_first = UINT32_MAX;
    uint32_t idx_second = UINT32_MAX;
    uint32_t sample_count = 0;
    float sum = 0;

    /// Checks if the current tile is a leaf by comparing if the two member indices are assigned.
    bool is_leaf() const;

    /// Determines if the leaf should be split by performing a One-Sample T-Test against the subdivision threshold.
    /// Important: is_leaf() should be called before to ensure this operation is only performed in a leaf!
    bool should_split(const int depth) const;
    
    /// Returns the split direction of this tile in a usable format.
    SplitDirection split_direction(const DepthCounter& depth_counter) const;

    /// Correctly updates the covariance, variance and mean deviation statistics.
    void update_statistics(const Sample& sample);

    /// Correctly updates sum and sample count. Must be called after update_statistics()!
    void update_sum(const Sample& sample);

    /// Returns the area-adjusted mean deviation of the current tile.
    float meandev(const int depth) const;

    /// Returns the area-adjusted variance of the current tile.
    float var(const int depth) const;

    /// Returns the covariance of the current tile given the splitting direction.
    float covar(const SplitDirection dir, const DepthCounter& depth_counter) const;

    /// Returns the absolute squared covariance used for determining the split direction.
    float adjusted_covar(const SplitDirection dir, const DepthCounter& depth_counter) const;

    /// Returns the mean of the luminance stored in this tile.
    float mean() const;

    /// Returns the ratio of the area of a tile (specified by the bounds) lying within the visible (i.e, sample-able) region.
    static float calc_visible_area_ratio(bool is_first, SplitDirection split_dir, Point2 x_bounds, Point2 y_bounds);
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
    std::vector<BinaryTile> tiles;

    Point2 x_bounds;
    Point2 y_bounds;
    float area_leaf_sum;

    /// Finds a tile based on the position stored in the uv parameter.
    BinaryTile& find_tile(const Point2& uv);

    /// Finds a tile based on the position stored in the uv parameter, but with the option to pass a DepthCounter to obtain the depth at which the tile is located.
    BinaryTile& find_tile(const Point2& uv, DepthCounter& counter);

    /// Stores a sample in a binary tiling.
    void insert(const Sample& sample);

    /// Calculate the area-adjusted sum of all means in leaf-nodes.
    float calc_leaf_sum();

    /// Warps a 2D coordinate in [0, 1]^2 to the range of the current tiling.
    Point2 warp_to_range(const Point2& uv) const;
};

/**
 * @brief Uppermost layer of the Binary Tile Coding (BTC).
 */
struct MTS_EXPORT_CORE BinaryTileCoding : public DataStructure {
    static float SUBDIV_THRESHOLD;
    static int MIN_SAMPLES;
    static int MAX_DEPTH;

    static Point2i tile_dims;

    ~BinaryTileCoding() { }

    void construct(DSArguments& init_data) override;
    void preprocess() override;
    void store(std::vector<Sample>& samples) override;
    void postprocess() override;

    Sample sample(Point2& pos) override;
    Float eval(Point2& pos) override;
    void wipe() override;
    DSType type() override;
    std::string name() override;
    int memory() override;

private:
    static RandomGen random;
    std::vector<BinaryTiling> tilings;
};

MTS_NAMESPACE_END 

#endif /* __DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_ */