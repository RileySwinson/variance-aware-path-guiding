#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_TILE_CODING_H_)
#define __DSCOMPARE_STRUCTURES_TILE_CODING_H_

#include <ds-compare/ds.h>
#include <Eigen/Core>

MTS_NAMESPACE_BEGIN

struct Tile {
    Float value = 0;
    int entries = 0;
};

typedef std::vector<Tile> Tiling;

struct MTS_EXPORT_CORE TileCoding : public DataStructure {
    ~TileCoding() { }

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
    std::vector<Tiling> tilings;
    std::vector<Float> guiding_map;

    int m_tiling_count = 4;
    Point2i m_tiling_dims; // x = width, y = height
    Sample::Mode m_mode;
    
    Float m_integral = 0;
    std::vector<Float> m_row_avgs;

    Float pdf(Point2& pos);
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_TILE_CODING_H_ */