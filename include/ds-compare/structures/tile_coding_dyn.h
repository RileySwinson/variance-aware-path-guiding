#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_TILE_CODING_DYNAMIC_H_)
#define __DSCOMPARE_STRUCTURES_TILE_CODING_DYNAMIC_H_

#include <ds-compare/ds.h>

MTS_NAMESPACE_BEGIN

struct DynamicTile {
    Float value = 0;
    std::array<uint32_t, 4> node_idx;
    std::vector<uint32_t> sample_idx;
};

typedef std::vector<DynamicTile> DynamicTiling;

struct MTS_EXPORT_CORE DynamicTileCoding : public DataStructure {
    ~DynamicTileCoding() { }

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
    std::vector<DynamicTiling> tilings;
    std::vector<Sample> samples;
    std::vector<Float> guiding_map;

    int m_tiling_count = 4;
    int m_tiling_dim = 16;
    int m_max_depth = 10;
    Float m_subdiv_thresh = 0.8;
    Sample::Mode m_mode;

    Float m_highest = 0;

    Float pdf(Point2& pos);
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_TILE_CODING_DYNAMIC_H_ */