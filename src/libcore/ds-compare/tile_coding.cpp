#include <ds-compare/structures/tile_coding.h>

MTS_NAMESPACE_BEGIN

RandomGen TileCoding::random = RandomGen();

Float Tile::area(int y, Point2i& inner)
{
    Float theta_step = M_PI / inner.y;
    Float d_phi = (2 * M_PI) / inner.x;
    Float d_theta = std::abs(std::cos(y * theta_step) - std::cos((y + 1) * theta_step));
    
    return (d_phi * d_theta);
}

void TileCoding::construct(DSArguments& init_data)
{
    SAssert(init_data.tc.tilings > 0 && init_data.tc.tiles_x > 0 && init_data.tc.tiles_y > 0);

    this->m_tiling_count = init_data.tc.tilings;
    this->m_tiling_dims = Point2i(init_data.tc.tiles_x, init_data.tc.tiles_y);
    this->m_mode = init_data.comparer.mode;
}

void TileCoding::preprocess()
{
    tilings.resize(this->m_tiling_count);
    for (auto& tiling : this->tilings)
    {
        tiling.resize(this->m_tiling_dims.x * this->m_tiling_dims.y);
    }
}

void TileCoding::store(std::vector<Sample>& samples)
{
    // Calculate offset
    Float tile_width = 1.0 / this->m_tiling_dims.x;
    Float tile_height = 1.0 / this->m_tiling_dims.y;

    Float overhead = 1.0 / this->m_tiling_count;
    if (this->m_tiling_count == 1) overhead = 0; // No offset shenanigans if we only have a single tile. Just span it over the whole thing.

    Point2 offset(tile_width * overhead, tile_height * overhead);

    // Precompute mapping ranges and factors
    std::vector<Point2> start_vals(this->m_tiling_count);
    std::vector<Float> x_slopes(this->m_tiling_count);
    std::vector<Float> y_slopes(this->m_tiling_count);
    for (int i = 0; i < this->m_tiling_count; ++i)
    {
        Point2 x_new(0 - (i * offset.x), 1 + ((this->m_tiling_count - 1 - i) * offset.x));
        Point2 y_new(0 - ((this->m_tiling_count - 1 - i) * offset.y), 1 + (i * offset.y));
    
        Float slope_x = 1.0 / (x_new.y - x_new.x);
        Float slope_y = 1.0 / (y_new.y - y_new.x);

        start_vals.at(i) = Point2(x_new.x, y_new.x);
        x_slopes.at(i) = slope_x;
        y_slopes.at(i) = slope_y;
    }

    // Store samples by mapping (sample range & offset) -> [0, 1) -> [0, dim(axis))
    for (const auto& sample : samples)
    {
        Point2 uv = Converter::spherical_to_uv(Point2f(sample.phi, sample.theta));

        for (int t_i = 0; t_i < this->m_tiling_count; ++t_i)
        {
            Tiling& tiling = this->tilings.at(t_i);

            Float warped_x = x_slopes.at(t_i) * (uv.x - start_vals.at(t_i).x);
            Float warped_y = y_slopes.at(t_i) * (uv.y - start_vals.at(t_i).y);

            Point2i index(
                warped_x * this->m_tiling_dims.x,
                warped_y * this->m_tiling_dims.y
            );

            int i = (index.y * this->m_tiling_dims.x) + index.x;
            Tile& tile = tiling.at(i);
            tile.sum += sample.value;
            tile.entries++;
        }
    }
}

void TileCoding::postprocess()
{
    /* Smush tilings to a single map */
    int x = this->m_tiling_dims.x * this->m_tiling_count;
    int y = this->m_tiling_dims.y * this->m_tiling_count;
    int total_overhead = this->m_tiling_count - 1;

    Point2i inner(x - total_overhead, y - total_overhead);
    int map_size = inner.x * inner.y;

    this->guiding_map.resize(map_size);

    // Track the total sum for normalization
    Float tile_area = 0;
    for (int i = 0; i < map_size; ++i)
    {
        Tile& t_final = this->guiding_map.at(i);

        const int pos_x = i % inner.x;
        const int pos_y = i / inner.x;

        // Calculate area to consider spherical trafo!
        if (pos_x == 0)
        {
            tile_area = Tile::area(pos_y, inner);
        }

        const Point2i base(0, total_overhead);
        for (int t_i = 0; t_i < this->tilings.size(); ++t_i)
        {
            Tiling& tiling = this->tilings.at(t_i);

            const Point2i pos(
                (base.x + t_i) + pos_x,
                (base.y - t_i) + pos_y
            );

            const int t_x = pos.x / this->m_tiling_count;
            const int t_y = pos.y / this->m_tiling_count;

            Tile& tile = tiling.at((t_y * this->m_tiling_dims.x) + t_x);
            if (tile.entries == 0) continue;

            t_final.sum += tile.sum;
            t_final.entries += tile.entries;
        }

        Float p_x = (t_final.sum / (t_final.entries * this->m_tiling_count));
        this->m_total_sum += p_x * tile_area;
    }

    /* Normalize & Precompute means */
    this->m_row_avgs.reserve(inner.y);

    float total_pdf_sum = 0.0f;
    for (int y = 0; y < inner.y; ++y)
    {
        float row_pdf_sum = 0.0f;
        for (int x = 0; x < inner.x; ++x)
        {
            int tile_i = (y * inner.x) + x;

            Tile& tile = this->guiding_map.at(tile_i);
            auto p_x = tile.sum / (tile.entries * this->m_total_sum);

            if (!(p_x > 0))
            {
                p_x = Epsilon;
            }

            row_pdf_sum += p_x;
        }

        total_pdf_sum += row_pdf_sum;
        this->m_row_avgs.push_back(row_pdf_sum / inner.x);
    }

    this->m_integral = total_pdf_sum / map_size;

    /* Clear tilings -- we only need the map */
    this->tilings.clear();
}

Sample TileCoding::sample(Point2& sample)
{
    int total_overhead = this->m_tiling_count - 1;
    int x_len = (this->m_tiling_dims.x * this->m_tiling_count) - total_overhead;
    int y_len = this->m_row_avgs.size();

    Float sum_y = 0;
    int y = 0;
    for (y = 0; y < y_len; ++y)
    {
        sum_y += this->m_row_avgs.at(y) / this->m_integral;
        if (sum_y / y_len >= sample.y) break;
    }
    if (y == y_len) y -= 1;

    Float sum_x = 0;
    int x = 0;
    for (x = 0; x < x_len; ++x)
    {
        int i = (y * x_len) + x;
        Tile& t = this->guiding_map.at(i);

        auto p_x = t.sum / (t.entries * this->m_total_sum);
        sum_x += p_x / this->m_row_avgs.at(y);

        if (sum_x / x_len >= sample.x) break;
    }
    if (x == x_len) x -= 1;

    Point2 tile_start(x / (Float) x_len, y / (Float) y_len);

    Point2 rng = random.next2D();
    Point2 uv(
        boost::algorithm::clamp(tile_start.x + rng.x * (1.0 / x_len), 0.0, 1.0 - Epsilon),
        boost::algorithm::clamp(tile_start.y + rng.y * (1.0 / y_len), 0.0, 1.0 - Epsilon)
    );

    Point2 spherical = Converter::uv_to_spherical(uv);

    Sample sample_data = {
        .value = 0,
        .pdf = pdf(uv),
        .theta = spherical.y,
        .phi = spherical.x
    };

    return sample_data;
}

Float TileCoding::eval(Point2& pos)
{
    return pdf(pos);
}

void TileCoding::wipe()
{
    this->guiding_map.clear();
    this->m_row_avgs.clear();
    this->m_integral = 0;
}

DSType TileCoding::type()
{
    return DSType::DS_TileCoding;
}

std::string TileCoding::name()
{
    return "Tile Coding";
}

Float TileCoding::pdf(Point2& pos)
{
    int x = this->m_tiling_dims.x * this->m_tiling_count;
    int y = this->m_tiling_dims.y * this->m_tiling_count;
    int total_overhead = this->m_tiling_count - 1;

    Point2i index(
        pos.x * (x - total_overhead),
        pos.y * (y - total_overhead)
    );

    Tile& tile = this->guiding_map.at((index.y * (x - total_overhead)) + index.x);
    auto p_x = tile.sum / (tile.entries * this->m_tiling_count * this->m_total_sum);

    return p_x;
}

int TileCoding::memory()
{
    // Calc size of tilings
    int tilings_size = 0;
    for (const auto& tiling : this->tilings)
    {
        // tile bytes * number of Tiles + vector bytes
        const int tiling_size = (sizeof(Tile) * tiling.capacity()) + sizeof(Tiling);
        tilings_size += tiling_size;
    }
    tilings_size += sizeof(this->tilings);

    // Calc size of map
    int map_size = sizeof(float) * this->guiding_map.capacity() + sizeof(this->guiding_map);

    return sizeof(this) + tilings_size + map_size;
}

MTS_NAMESPACE_END