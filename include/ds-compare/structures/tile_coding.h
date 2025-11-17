#pragma once
#if !defined(__DSCOMPARE_STRUCTURES_TILE_CODING_H_)
#define __DSCOMPARE_STRUCTURES_TILE_CODING_H_

#include <ds-compare/ds.h>

MTS_NAMESPACE_BEGIN

struct Tile {
    float sum = 0;

    static Float area(const Point2& y_pos_norm, const int map_tiles_x);
};

typedef std::vector<Tile> Tiling;

struct GuidingMap {
    Tiling tiles;
    Point2i dims;

    /// Returns the weighted sum stored at the passed in position. If the position contains no value, Epsilon is returned instead.
    float get(int pos);
};

struct MTS_EXPORT_CORE TileCoding : public DataStructure {
    ~TileCoding() { }

    void construct_impl(DSArguments& init_data) override;
    void preprocess_impl() override;
    void store_impl(Sample& sample) override;
    void postprocess_impl(bool last_iteration = false) override;

    Sample sample_impl(Point2& pos) override;
    Float eval_impl(Point2& pos) override;
    void wipe_impl() override;
    DSType type_impl() override;
    std::string name_impl() override;
    int memory_impl() override;

private:
    std::vector<Tiling> tilings;
    GuidingMap guiding_map;

    int m_tiling_count = 4;
    Point2i m_tiling_dims; // x = width, y = height
    Sample::Mode m_mode;
    TCParams::Transformation m_transform_mode;

    std::vector<Point2> m_tiling_origin;
    Point2 m_tiling_len;
    
    Float m_integral = 0;
    std::vector<Float> m_cdf;
    Float m_total_sum = 0;

    RandomGen random;

    Float pdf(Point2& pos);
    void build_map();
    void calc_cdf();
    void empty_tilings();

    template <typename T>
    T transform(T value, bool inverse = false, const std::function<T(T)>& f = nullptr);
};

MTS_NAMESPACE_END

#endif /* __DSCOMPARE_STRUCTURES_TILE_CODING_H_ */