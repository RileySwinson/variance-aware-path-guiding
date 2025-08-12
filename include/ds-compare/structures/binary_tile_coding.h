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

enum TileType {
    Node,
    Leaf
};

/**
 * @brief Helper struct to track tile data.
 */
struct TileTracker {
    int depth = 0;
    Point2 x_bounds;
    Point2 y_bounds;

    inline void increment()
    {
        this->depth++;
    }

    inline int depth() const
    {
        return this->depth;
    }

    inline void boundaries(const Point2 x, const Point2 y)
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

/**
 * @brief Tile in a tiling.
 * 
 * The most low-level entity in the 'Binary Tile Coding' data structure, storing radiance information
 * in a specific area of the sample space (and beyond). Can be both leaf and non-leaf (node), determined
 * by the should_split() method that returns true if the collected sample statistics implicate that the
 * tile should be split. To access a child, use the tiles vector in combination with the child indices.
 */
struct BinaryTile {
    struct TileData {
        uint32_t sample_count : 31;
        TileType tile_type    : 1;

        TileData() : sample_count(0), tile_type(Leaf) { };
    };

    struct NodeData {
        /* ==== Traversal Data (20 bytes) ==== */
        std::array<float, 2> power = { 0.0f, 0.0f };
        std::array<uint32_t, 2> children = { UINT32_MAX, UINT32_MAX };
        Direction split_direction = Horizontal;
    };

    struct LeafData {
        /* ==== Statistics (24 bytes) ==== */
        Point2f cov = Point2f(0.0f);
        Point2f sample_mean = Point2f(0.0f);
        float m2 = 0;
        float diff_sum = 0;
    };

    TileData data;
    float sum = 0.0f;
    union {
        NodeData node;
        LeafData leaf;
    };

    BinaryTile(const TileType type)
    {
        this->data = TileData();
        if (type == TileType::Node)
            this->node = NodeData();
        else
            this->leaf = LeafData();
    };

    /// Returns whether the current tile is a leaf.
    bool is_leaf() const;

    /// Determines if a leaf tile should be split by performing a One-Sample T-Test against the subdivision threshold.
    bool should_split(const TileTracker& tracker) const;
    
    /// Returns the split direction of a leaf tile by calculating the absolute covariance in both x and y direction.
    Direction split_direction(const TileTracker& tracker) const;

    /// Correctly updates the covariance, variance and mean deviation statistics.
    void update_statistics(const Sample& sample);

    /// Correctly updates sum and sample count. Must be called after update_statistics()!
    void update_sum(const Sample& sample);

    /// Returns the area-adjusted mean deviation of the current leaf.
    float meandev(const TileTracker& tracker) const;

    /// Returns the area-adjusted variance of the current leaf.
    float var(const TileTracker& tracker) const;

    /// Returns the mean of the radiance stored in this tile.
    float mean() const;

    /// Returns the normalized area of this tile.
    float area(const TileTracker& tracker) const;

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
    std::vector<BinaryTile> tiles;
    Point2 x_bounds;
    Point2 y_bounds;

    float total_base_mean = 0.0f;
    std::vector<float> row_means;
    std::vector<float> tile_areas;

    /// Finds a tile at a given position within the bounds of this tiling. If find_empty is false, the parent will be returned in case both children are empty.
    BinaryTile& find_tile(const Point2& pos, bool find_empty = true);

    /// Finds a tile at a given position within the bounds of this tiling, but with the option to pass a TileTracker to obtain further data about the tile. If find_empty is false, the parent will be returned in case both children are empty.
    BinaryTile& find_tile(const Point2& pos, TileTracker& counter, bool find_empty = true);

    /// Stores a sample in a binary tiling.
    void insert(const Sample& sample);

    /// Utility function to recursively add the power (mean * area) from all children to their parent tiles.
    float recurse_statistics(BinaryTile& curr_tile, TileTracker tracker, float& leaf_sum);

    /// Utility function to calculate the planar (!) boundaries of a base tile.
    std::pair<Point2, Point2> base_tile_bounds(int x, int y) const;

    /// Utility function to obtain the x and y position of a base tile by sampling from a CDF.
    Point2i base_tile_pos(const Point2& pos);

    /// Checks if the children of the passed in tile are both empty, i.e., their mean equals 0.
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
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_ */