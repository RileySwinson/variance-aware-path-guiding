#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_)
#define __DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_

#include <ds-compare/ds.h>

MTS_NAMESPACE_BEGIN

// Set this to true to enable multithreading which should speed things up a fair bit IFF a somewhat significant number of samples are proccessed at once.
#define BTC_MULTITHREADING false

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
    uint32_t index = 0;
    int depth = 0;
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

    inline void set_index(const uint32_t index)
    {
        this->index = index;
    }

    inline Point2 clamped(const Direction dir) const
    {
        Point2 bounds = (dir == Horizontal) ? this->x_bounds : this->y_bounds;
        return Point2(
            boost::algorithm::clamp(bounds.x, 0.0, 1.0),
            boost::algorithm::clamp(bounds.y, 0.0, 1.0)
        );
    }
};

/**
 * @brief Tracks the statistics of a tile.
 *
 * Similarly to how each tile is stored, a tiling also stores a vector of statistics objects, one per non-leaf tile.
 * Upon sample storage, the respective stats object is fetched and its statistics are updated via a 1-pass
 * Welfold algorithm.
 * 
 * The stats primarily decide two things:
 * 1) Checking if a tile should be split, by conducting a Student's 1-sample t-test
 * 2) Checking how a tile should be split, by looking at the covariances in both directions (phi, theta)
 * 
 * Keep in mind that after the final learning iteration (i.e., when the learning phase is over) the stats vector
 * gets cleared by a user-initialized postprocess() call.
 */
struct BTStatistics {
    /// Correctly updates all required statistics.
    void update(const Sample& sample, const float tile_sum);

    /// Returns the split direction of a leaf tile by calculating the absolute covariance in both x (phi) and y (theta) direction.
    Direction split_direction(const BTTracker& tracker) const;

    /// Determines if a leaf tile should be split by performing a one-sample t-test, trying to see if we can reject our null hypothesis that the MD is significantly different from 0.
    bool should_split() const;

private:
    // The covariance in x and y direction to determine how to split.
    Point2f cov = Point2f(0.0f);
    // Sample mean required for covariance updates.
    Point2f sample_mean = Point2f(0.0f);
    // Sum of squared differences of the absolute deviations from their mean (used for variance).
    float m2 = 0.0f;
    // Sum of absolute differences used for MD metric.
    float diff_sum = 0.0f;
    // Number of samples that arrived in a tile.
    uint32_t sample_count = 0;

    /// Returns the mean deviation of the current leaf.
    float meandev() const;

    /// Returns the variance *of the mean deviations* of the current leaf.
    float var() const;
};

/**
 * @brief Tile in a tiling.
 * 
 * The most low-level entity in the 'Binary Tile Coding' data structure, storing radiance information
 * in a specific area of the sample space (and beyond). Can be both leaf and non-leaf (node), determined
 * by the is_leaf() method that returns true if both children are set to valid indices.
 * To access a child, use the tiles vector in combination with the child indices.
 */
struct BinaryTile {
    float sum = 0.0f;
    std::array<float, 2> power = { 0.0f, 0.0f };
    std::array<uint32_t, 2> children = { UINT32_MAX, UINT32_MAX };
    Direction split_direction = Horizontal;

    /// Returns whether the current tile is a leaf.
    bool is_leaf() const;

    /// Returns the normalized area of this tile.
    float area(const BTTracker& tracker);

    /// Update the weighted sum stored in this tile.
    void update(const Sample& sample);
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
    std::vector<BTStatistics> stats;
    Point2 x_bounds;
    Point2 y_bounds;

    /* Information used for marginal and conditional sampling (CDF-based sampling) */
    float total_power = 0.0f;
    std::vector<float> row_means;
    std::vector<float> tile_areas;

    /// Finds a tile at a given position within the bounds of this tiling. If find_empty is false, the parent will be returned in case both children are empty.
    BinaryTile& find_tile(const Point2& pos, bool find_empty = true);

    /// Finds a tile at a given position within the bounds of this tiling, but with the option to pass a tracker to obtain further data about the tile. If find_empty is false, the parent will be returned in case both children are empty.
    BinaryTile& find_tile(const Point2& pos, BTTracker& counter, bool find_empty = true);

    /// Stores a sample in a binary tiling.
    void insert(const Sample& sample);

    /// Utility function to recursively add the power (weighted sum * area) from all children to their parent tiles.
    float recurse_statistics(BinaryTile& curr_tile, BTTracker tracker, float& leaf_sum);

    /// Utility function to calculate the planar (!) boundaries of a base tile.
    std::pair<Point2, Point2> base_tile_bounds(int x, int y) const;

    /// Utility function to obtain the 2D index of a base tile.
    Point2i base_tile_pos(const Point2& pos);

    /// Checks if the children of the passed in tile are both empty, i.e., their weighted sum equals 0.
    inline bool children_empty(BinaryTile& tile) const;

    /// Transforms a value to the right domain (e.g., spherical, cosine). The user may specify a custom transformation by providing a lambda as last parameter.
    template <typename T>
    T transform(T value, bool inverse = false, const std::function<T(T)>& f = nullptr);
};

/**
 * @brief Uppermost layer of the Binary Tile Coding (BTC).
 */
struct MTS_EXPORT_CORE BinaryTileCoding : public DataStructure {
    static TTable::CI CI;
    static TCParams::Transformation TRANSFORM_MODE;
    static int MAX_DEPTH;
    static Point2i TILE_DIMS;

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

    void build();
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_ */