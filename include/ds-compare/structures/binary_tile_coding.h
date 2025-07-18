#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_)
#define __DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_

#include <ds-compare/ds.h>

MTS_NAMESPACE_BEGIN

/**
 * @brief Used to specify the direction of an entity, e.g. a split, amongst other things.
 * 
 * In the case of a split:
 * HORIZONTAL = divides horizontally, i.e. the split goes from top to bottom.
 * VERTICAL = divides vertically, i.e. the split goes from left to right.
 */
enum Direction {
    HORIZONTAL,
    VERTICAL
};

/**
 * @brief Used to specify how a tile is split.
 * 
 * PLANAR = the tile will be split as if on a 2D plane.
 * SPHERICAL = the tile will be split with spherical trafo in mind.
 */
enum SplitBehavior {
    PLANAR,
    SPHERICAL
};

/**
 * @brief Helper struct to track tile data.
 */
struct TileTracker {
    Point2i splits = Point2i(0);
    Direction last = HORIZONTAL;
    Point2 x_bounds;
    Point2 y_bounds;
    bool before_split = true;

    inline int depth() const
    {
        return this->splits.x + this->splits.y;
    }

    inline void side(const bool decision)
    {
        this->before_split = decision;
    }

    inline void boundaries(const Point2 x, const Point2 y)
    {
        this->x_bounds = x;
        this->y_bounds = y;
    }

    inline void increment(const Direction split)
    {
        ((split == HORIZONTAL) ? this->splits.x : this->splits.y)++;
        this->last = split;
    }

    inline Point2 clamped(const Direction dir) const
    {
        Point2 bounds = (dir == HORIZONTAL) ? this->x_bounds : this->y_bounds;
        return Point2(
            boost::algorithm::clamp(bounds.x, 0.0, 1.0),
            boost::algorithm::clamp(bounds.y, 0.0, 1.0)
        );
    }
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
    bool should_split(const TileTracker& tracker) const;
    
    /// Returns the split direction of this tile by calculating the absolute covariance in both x and y direction.
    Direction split_direction(const TileTracker& tracker) const;

    /// Correctly updates the covariance, variance and mean deviation statistics.
    void update_statistics(const Sample& sample);

    /// Correctly updates sum and sample count. Must be called after update_statistics()!
    void update_sum(const Sample& sample);

    /// Returns the area-adjusted mean deviation of the current tile.
    float meandev(const TileTracker& tracker) const;

    /// Returns the area-adjusted variance of the current tile.
    float var(const TileTracker& tracker) const;

    /// Returns the mean of the luminance stored in this tile.
    float mean() const;

    /// Returns the normalized area of this tile.
    float area(const TileTracker& tracker) const;
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

    /// Finds a tile at a given position within the bounds of this tiling. If find_empty is false, the parent will be returned in case both children are empty.
    BinaryTile& find_tile(const Point2& pos, bool find_empty = true);

    /// Finds a tile at a given position within the bounds of this tiling, but with the option to pass a TileTracker to obtain further data about the tile. If find_empty is false, the parent will be returned in case both children are empty.
    BinaryTile& find_tile(const Point2& pos, TileTracker& counter, bool find_empty = true);

    /// Stores a sample in a binary tiling.
    void insert(const Sample& sample);

    /// Utility function to recursively add the sum and sample count statistics from child tiles to parent tiles.
    std::pair<uint32_t, float> recurse_statistics(BinaryTile& curr_tile, TileTracker tracker, float& leaf_sum);

    /// Utility function to obtain the x and y position of a base tile by sampling from a CDF.
    Point2i base_tile_pos_from_cdf(const Point2& pos) const;

    /// Calculates the split position of a tile depending on the used split behavior.
    Float split_position(const Point2& bounds, const Direction split_dir) const;

    /// Checks if the children of the passed in tile are both empty, i.e., their mean equals 0.
    inline bool children_empty(BinaryTile& tile) const;
};

/**
 * @brief Uppermost layer of the Binary Tile Coding (BTC).
 */
struct MTS_EXPORT_CORE BinaryTileCoding : public DataStructure {
    static int MAX_DEPTH;
    static TTable::CI ci;
    static SplitBehavior split_behavior;
    static Point2i tile_dims;
    float leaf_sum = 0.0f;

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