#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_)
#define __DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_

#include <ds-compare/ds.h>

MTS_NAMESPACE_BEGIN

enum SplitDirection {
    HORIZONTAL,
    VERTICAL
};

struct Bitfield32 {
    unsigned int split : 1; // how a tile is split -- 0 = horizontal, 1 = vertical
    unsigned int sample_count : 31; // number of samples in a tile
};

struct BinaryTile {
    Point2f cov;
    Point2f sample_mean;
    float value = 0;

    uint32_t idx_first = -1;
    uint32_t idx_second = -1;

    Bitfield32 meta_data;

    bool is_leaf() const;
    bool should_split() const;
    SplitDirection split_direction() const;
    void update_covariance(Sample& sample);
    void update_value(Sample& sample);
    float covar(SplitDirection dir) const;
};

struct BinaryTiling {
    std::vector<BinaryTile> tiles;

    Point2f start_vals;
    Point2f factors;

    BinaryTile* find_tile(Point2& uv, Point2i& tile_dims, int& depth);
    void insert(Sample& sample, Point2i& tile_dims);
};

struct MTS_EXPORT_CORE BinaryTileCoding : public DataStructure {
    static constexpr float SUBDIV_THRESHOLD = 0.10f;
    static constexpr int MIN_SAMPLES = 100;
    static constexpr int MAX_DEPTH = 10;

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
    std::vector<BinaryTiling> tilings;
    Point2i tile_dims;
};

MTS_NAMESPACE_END 

#endif /* __DSCOMPARE_STRUCTURES_BINARY_TILE_CODING_H_ */