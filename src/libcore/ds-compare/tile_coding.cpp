#include <ds-compare/structures/tile_coding.h>

MTS_NAMESPACE_BEGIN

Float Tile::area(int y, Point2i& inner)
{
    Float theta_step = M_PI / inner.y;
    Float d_phi = (2 * M_PI) / inner.x;
    Float d_theta = std::abs(std::cos(y * theta_step) - std::cos((y + 1) * theta_step));
    
    return (d_phi * d_theta);
}

float GuidingMap::mean(int pos)
{
    Tile& t = this->tiles[pos];

    if (t.sum == 0 || t.entries == 0)
    {
        return Epsilon;
    }

    return t.sum / t.entries;
}

void TileCoding::construct(DSArguments& init_data)
{
    SAssert(init_data.tc.tilings > 0 && init_data.tc.tiles_x > 0 && init_data.tc.tiles_y > 0);

    this->m_tiling_dims = Point2i(init_data.tc.tiles_x, init_data.tc.tiles_y);
    this->m_mode = init_data.comparer.mode;
    this->m_tiling_count = init_data.tc.tilings;
}

void TileCoding::preprocess()
{
    reset_tilings();

    Point2i inner(
        (this->m_tiling_dims.x * this->m_tiling_count) - (this->m_tiling_count - 1),
        (this->m_tiling_dims.y * this->m_tiling_count) - (this->m_tiling_count - 1)
    );

    this->guiding_map.tiles.resize(inner.x * inner.y);
    this->guiding_map.dims = inner;
}

void TileCoding::store(std::vector<Sample>& samples)
{
    const auto dims_x = this->m_tiling_dims.x;
    const auto dims_y = this->m_tiling_dims.y;

    // Calculate offset
    Float tile_width = 1.0 / dims_x;
    Float tile_height = 1.0 / dims_y;

    Float overhead = 1.0 / this->m_tiling_count;
    if (this->m_tiling_count == 1) overhead = 0; // No offset shenanigans if we only have a single tile. Just span it over the whole thing.

    Point2 shift(tile_width * overhead, tile_height * overhead);
    Float x_len = 1 + ((this->m_tiling_count - 1) * shift.x);
    Float y_len = 1 + ((this->m_tiling_count - 1) * shift.y);

    // Store samples by mapping (sample range & offset) -> [0, 1) -> [0, dim(axis))
    for (const auto& sample : samples)
    {
        Point2 uv = Converter::spherical_to_uv(Point2(sample.phi, sample.theta));

        for (int ti = 0; ti < this->m_tiling_count; ++ti)
        {
            Tiling& tiling = this->tilings[ti];

            Point2 t_origin(0 - (ti * shift.x), 0 - (ti * shift.y));

            Float x_warped = Converter::map(uv.x).from({ t_origin.x, t_origin.x + x_len }).to({ 0.0, 1.0 });
            Float y_warped = Converter::map(uv.y).from({ t_origin.y, t_origin.y + y_len }).to({ 0.0, 1.0 });

            Point2i index(
                x_warped * dims_x,
                y_warped * dims_y
            );

            if (index.x == dims_x) index.x--;
            if (index.y == dims_y) index.y--;

            int i = (index.y * dims_x) + index.x;
            Tile& tile = tiling[i];
            tile.sum += sample.value / sample.pdf;
            tile.entries++;
        }
    }

    build_map();
    calc_cdf();
    reset_tilings();
}

void TileCoding::postprocess()
{
    /* Wipe tilings completely -- we only need the map */
    this->tilings.clear();
    this->tilings.shrink_to_fit();
}

Sample TileCoding::sample(Point2& sample)
{
    int total_overhead = this->m_tiling_count - 1;
    int x_len = (this->m_tiling_dims.x * this->m_tiling_count) - total_overhead;
    int y_len = this->m_row_avgs.size();

    Float sum_y = 0; int y = 0;
    for (y = 0; y < y_len; ++y)
    {
        sum_y += this->m_row_avgs[y];
        if (sum_y / (y_len * this->m_integral) >= sample.y) break;
    }
    if (y == y_len) y -= 1;

    Float sum_x = 0; int x = 0;
    for (x = 0; x < x_len; ++x)
    {
        int i = (y * x_len) + x;
        sum_x += this->guiding_map.mean(i);
        if (sum_x / (x_len * this->m_total_sum * this->m_row_avgs[y]) >= sample.x) break;
    }
    if (x == x_len) x -= 1;

    Point2 rng = this->random.next2D();

    Point2 x_bounds(x / (Float) x_len, (x + 1) / (Float) x_len);
    Point2 y_bounds(y / (Float) y_len, (y + 1) / (Float) y_len);
    y_bounds.x = std::cos(y_bounds.x * M_PI);
    y_bounds.y = std::cos(y_bounds.y * M_PI);

    Point2 uv(
        x_bounds.x + rng.x * (x_bounds.y - x_bounds.x),
        y_bounds.x + rng.y * (y_bounds.y - y_bounds.x)
    );
    uv.y = std::acos(uv.y) * INV_PI;

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
    this->guiding_map.tiles.clear();
    this->m_row_avgs.clear();
    this->m_integral = 0;
    this->m_total_sum = 0;
}

DSType TileCoding::type()
{
    return DSType::DS_TileCoding;
}

std::string TileCoding::name()
{
    return "Tile Coding";
}

int TileCoding::memory()
{
    size_t size_self = sizeof(*this);
    size_t size_tilings = this->tilings.capacity() * sizeof(Tiling);

    size_t size_tiles = 0;
    for (const auto& tiling : this->tilings)
    {
        size_tiles += tiling.capacity() * sizeof(Tile);
    }

    size_t size_map = this->guiding_map.tiles.capacity() * sizeof(Tile);

    return size_self + size_tilings + size_tiles + size_map;
}

Float TileCoding::pdf(Point2& pos)
{
    int x = this->m_tiling_dims.x * this->m_tiling_count;
    int y = this->m_tiling_dims.y * this->m_tiling_count;
    int x_tiles = x - (this->m_tiling_count - 1);
    int y_tiles = y - (this->m_tiling_count - 1);

    Point2i index(pos.x * x_tiles, pos.y * y_tiles);
    if (index.x == x_tiles) index.x--;
    if (index.y == y_tiles) index.y--;

    int i = (index.y * x_tiles) + index.x;
    float t_mu = this->guiding_map.mean(i);
    return t_mu / (this->m_tiling_count * this->m_total_sum);
}

void TileCoding::build_map()
{
    this->m_total_sum = 0;

    Tiling& map = this->guiding_map.tiles;
    Point2i dims = this->guiding_map.dims;
    const int map_size = dims.x * dims.y;

    // Track the total sum for normalization
    Float tile_area = 0;
    for (int i = 0; i < map_size; ++i)
    {
        const int pos_x = i % dims.x;
        const int pos_y = i / dims.x;

        // Calculate area to consider spherical trafo!
        if (pos_x == 0)
        {
            tile_area = Tile::area(pos_y, dims);
        }

        const Point2i base(0, this->tilings.size() - 1);
        for (size_t t_i = 0; t_i < this->tilings.size(); ++t_i)
        {
            Tiling& tiling = this->tilings[t_i];

            const int t_x = (t_i + pos_x) / this->tilings.size();
            const int t_y = (t_i + pos_y) / this->tilings.size();

            Tile& tile = tiling[(t_y * this->m_tiling_dims.x) + t_x];
            if (tile.entries == 0) continue;

            map[i].sum += tile.sum;
            map[i].entries += tile.entries;
        }

        Float p_x = this->guiding_map.mean(i) / this->tilings.size();
        this->m_total_sum += p_x * tile_area;
    }
}

void TileCoding::calc_cdf()
{
    const Point2i dims = this->guiding_map.dims;

    /* Normalize & Precompute means */
    this->m_row_avgs = std::vector<Float>(dims.y);
    this->m_integral = 0.0;

    float total_pdf_sum = 0.0f;
    for (int y = 0; y < dims.y; ++y)
    {
        float row_pdf_sum = 0.0f;
        for (int x = 0; x < dims.x; ++x)
        {
            int tile_i = (y * dims.x) + x;
            auto p_x = this->guiding_map.mean(tile_i) / this->m_total_sum;
            row_pdf_sum += p_x;
        }

        total_pdf_sum += row_pdf_sum;
        this->m_row_avgs[y] = row_pdf_sum / dims.x;
    }

    this->m_integral = total_pdf_sum / (dims.x * dims.y);
}

void TileCoding::reset_tilings()
{
    this->tilings = std::vector<Tiling>(this->m_tiling_count);
    for (auto& tiling : this->tilings)
    {
        tiling.resize(this->m_tiling_dims.x * this->m_tiling_dims.y);
    }
}

MTS_NAMESPACE_END