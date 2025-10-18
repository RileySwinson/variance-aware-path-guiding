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

    inline Point2 clamped(const Direction dir) const
    {
        Point2 bounds = (dir == Horizontal) ? this->x_bounds : this->y_bounds;
        return Point2(
            boost::algorithm::clamp(bounds.x, 0.0, 1.0),
            boost::algorithm::clamp(bounds.y, 0.0, 1.0)
        );
    }
};

enum TileType {
    Internal,
    Leaf
};

struct BTBitfield {
    uint32_t sample_count : 31;
    TileType tile_type : 1;

    BTBitfield() : sample_count(0), tile_type(Leaf) {};
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
    Point2f cov = Point2f(0.0f);
    // Sample mean required for covariance updates.
    Point2f sample_mean = Point2f(0.0f);
    // Sum of squared distances from the mean.
    float m2 = 0.0f;
    // Sum of absolute deviations used for MD metric.
    float ad_sum = 0.0f;
};

/**
 * @brief Stores relevant data about child tiles.
 */
struct InternalNode {
    // The total radiance in each child (including sub-children).
    std::array<float, 2> power = { 0.0f, 0.0f };
    // The indices of the children. Use in combination with the tiles vector.
    std::array<uint32_t, 2> children = { UINT32_MAX, UINT32_MAX };
    // How the children are positioned.
    Direction split_direction = Horizontal;
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
struct BinaryTile {
    // Total radiance that arrived in this tile.
    float sum = 0.0f;
    // Sample count + tile flag
    BTBitfield data;

    union {
        InternalNode internal;
        LeafNode leaf;
    };

    BinaryTile(const TileType type)
    {
        this->data = BTBitfield();
        if (type == TileType::Internal)
        {
            this->internal = InternalNode();
        }
        else
        {
            this->leaf = LeafNode();
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

    /// Returns the normalized area of this tile.
    float area(const BTTracker& tracker);

    /// Returns whether the current tile is a leaf.
    inline bool is_leaf() const;

private:
    BinaryTile() = default;
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

    /* Information used for marginal and conditional sampling (CDF-based sampling) */
    float total_power = 0.0f;
    std::vector<float> row_means;
    std::vector<float> tile_areas;

    /// Finds a tile at a given position within the bounds of this tiling. If find_empty is false, the parent will be returned in case both children are empty.
    BinaryTile& find_tile(const Point2& pos, bool find_empty = true);

    /// Finds a tile at a given position within the bounds of this tiling, but with the option to pass a tracker to obtain further data about the tile. If find_empty is false, the parent will be returned in case both children are empty.
    BinaryTile& find_tile(const Point2& pos, BTTracker& counter, bool find_empty = true);

    /// Stores a sample in a binary tiling.
    void insert(const Sample& sample, uint32_t& split_count);

    /// Utility function to recursively add the power (weighted sum * area) from all children to their parent tiles.
    float recurse_statistics(BinaryTile& curr_tile, BTTracker tracker, float& leaf_sum);

    /// Utility function to calculate the planar (!) boundaries of a base tile.
    inline std::pair<Point2, Point2> base_tile_bounds(int x, int y) const;

    /// Utility function to obtain the 2D index of a base tile.
    inline Point2i base_tile_pos(const Point2& pos);

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
    static TTable::CI                   CI;
    static TCParams::Transformation     TRANSFORM_MODE;
    static uint32_t                     MAX_DEPTH;
    static uint32_t                     MAX_SPLITS;
    static Point2i                      TILE_DIMS;
    static float                        TILING_EXCESS;

    float leaf_sum = 0.0f;
    uint32_t split_count = 0;

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