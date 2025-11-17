#include <ds-compare/structures/tile_coding.h>

MTS_NAMESPACE_BEGIN

Float Tile::area(const Point2& y_pos_norm, const int map_tiles_x)
{
    Float d_theta = std::cos(M_PI * y_pos_norm.x) - std::cos(M_PI * y_pos_norm.y);
    Float d_phi = (2 * M_PI) / map_tiles_x;
    
    return d_theta * d_phi;
}

float GuidingMap::get(int pos)
{
    Tile& t = this->tiles[pos];

    if (t.sum == 0)
    {
        return Epsilon;
    }

    return t.sum;
}

void TileCoding::construct_impl(DSArguments& init_data)
{
    SAssert(init_data.tc.tilings > 0 && init_data.tc.tiles_x > 0 && init_data.tc.tiles_y > 0);

    this->m_tiling_dims = Point2i(init_data.tc.tiles_x, init_data.tc.tiles_y);
    this->m_mode = init_data.comparer.mode;
    this->m_tiling_count = init_data.tc.tilings;
    this->m_transform_mode = init_data.tc.transformation_mode;
}

void TileCoding::preprocess_impl()
{
    empty_tilings();

    Point2i inner(
        (this->m_tiling_dims.x * this->m_tiling_count) - (this->m_tiling_count - 1),
        (this->m_tiling_dims.y * this->m_tiling_count) - (this->m_tiling_count - 1)
    );

    this->guiding_map.tiles.resize(inner.x * inner.y);
    this->guiding_map.dims = inner;

    // Calculate offset
    int total_shifts = this->m_tiling_count - 1;
    Point2 shift(
        1.0 / (this->m_tiling_dims.x * this->m_tiling_count - total_shifts),
        1.0 / (this->m_tiling_dims.y * this->m_tiling_count - total_shifts)
    );

    this->m_tiling_len = Point2(
        1 + (total_shifts * shift.x),
        1 + (total_shifts * shift.y)
    );

	this->m_tiling_origin.reserve(this->m_tiling_count);
    for (int ti = 0; ti < this->m_tiling_count; ++ti)
    {
        Point2 t_origin(0 - (ti * shift.x), 0 - (ti * shift.y));
        this->m_tiling_origin.push_back(t_origin);
    }
}

void TileCoding::store_impl(Sample& sample)
{
    const auto dims_x = this->m_tiling_dims.x;
    const auto dims_y = this->m_tiling_dims.y;

    Point2 uv = Converter::spherical_to_uv(Point2(sample.phi, sample.theta));
    uv.y = transform(uv.y, true);

    for (int ti = 0; ti < this->m_tiling_count; ++ti)
    {
		Point2& t_origin = this->m_tiling_origin[ti];
        Float x_warped = Converter::map(uv.x).from({ t_origin.x, t_origin.x + this->m_tiling_len.x }).to({ 0.0, 1.0 });
        Float y_warped = Converter::map(uv.y).from({ t_origin.y, t_origin.y + this->m_tiling_len.y }).to({ 0.0, 1.0 });

        Point2i index(
            x_warped * dims_x,
            y_warped * dims_y
        );

        if (index.x == dims_x) index.x--;
        if (index.y == dims_y) index.y--;

        Tiling& tiling = this->tilings[ti];
        int i = (index.y * dims_x) + index.x;
        Tile& tile = tiling[i];
        tile.sum += sample.value / sample.pdf;
	}
}

void TileCoding::postprocess_impl(bool last_iteration)
{
	build_map();
	calc_cdf();
    empty_tilings();

    /* Wipe tilings completely -- we only need the map */
    if (last_iteration)
    {
        this->tilings.clear();
        this->tilings.shrink_to_fit();
    }
}

Sample TileCoding::sample_impl(Point2& sample)
{
    if (this->m_total_sum <= 0)
    {
        Point2 rng = this->random.next2D();

        Sample empty;
        empty.value = 0;
        empty.pdf = 1.0 / (4 * M_PI);
        empty.phi = 2 * M_PI * rng.x;
        empty.theta = std::acos(1 - 2 * rng.y);

		return empty;
    }

    int total_shifts = this->m_tiling_count - 1;
    int x_len = (this->m_tiling_dims.x * this->m_tiling_count) - total_shifts;
    int y_len = this->m_cdf.size();

    Float sum_y = 0; int y = 0;
    for (y = 0; y < y_len; ++y)
    {
        sum_y += this->m_cdf[y];
        if (sum_y / (y_len * this->m_integral) >= sample.y) break;
    }
    if (y == y_len) y -= 1;

    Float sum_x = 0; int x = 0;
    for (x = 0; x < x_len; ++x)
    {
        int i = (y * x_len) + x;
        sum_x += this->guiding_map.get(i);
        if (sum_x / (x_len * this->m_total_sum * this->m_cdf[y]) >= sample.x) break;
    }
    if (x == x_len) x -= 1;

    Point2 rng = this->random.next2D();

    Point2 x_bounds(x / (Float) x_len, (x + 1) / (Float) x_len);
    Point2 y_bounds(y / (Float) y_len, (y + 1) / (Float) y_len);
    y_bounds.x = std::cos(transform(y_bounds.x) * M_PI);
    y_bounds.y = std::cos(transform(y_bounds.y) * M_PI);

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

Float TileCoding::eval_impl(Point2& pos)
{
    return pdf(pos);
}

void TileCoding::wipe_impl()
{
    this->guiding_map.tiles.clear();
    this->m_cdf.clear();
    this->m_integral = 0;
    this->m_total_sum = 0;
}

DSType TileCoding::type_impl()
{
    return DSType::DS_TileCoding;
}

std::string TileCoding::name_impl()
{
    return "Tile Coding";
}

int TileCoding::memory_impl()
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

    Point2 warped_pos(pos.x, transform(pos.y, true));

    Point2i index(
        warped_pos.x * x_tiles,
        warped_pos.y * y_tiles
    );

    if (index.x == x_tiles) index.x--;
    if (index.y == y_tiles) index.y--;

    int i = (index.y * x_tiles) + index.x;
    float t_mu = this->guiding_map.get(i);
    return t_mu / this->m_total_sum;
}

void TileCoding::build_map()
{
    this->m_total_sum = 0;

    Tiling& map = this->guiding_map.tiles;
    Point2i dims = this->guiding_map.dims;
    const int map_size = dims.x * dims.y;

    // Track the total sum for normalization
    Float tile_area = 0.0;
    for (int i = 0; i < map_size; ++i)
    {
        const int pos_x = i % dims.x;
        const int pos_y = i / dims.x;

        // Calculate area to consider spherical trafo!
        if (pos_x == 0)
        {
            if (this->m_transform_mode == TCParams::Transformation::Spherical)
            {
                tile_area = (4 * M_PI) / map_size;
            }
            else
            {
                Point2 y_bounds_norm(
                    transform(pos_y / (Float) dims.y),
                    transform((pos_y + 1) / (Float) dims.y)
                );

                tile_area = Tile::area(y_bounds_norm, dims.x);
            }
        }

        for (size_t t_i = 0; t_i < this->tilings.size(); ++t_i)
        {
            Tiling& tiling = this->tilings[t_i];

            const int t_x = (t_i + pos_x) / this->tilings.size();
            const int t_y = (t_i + pos_y) / this->tilings.size();

            Tile& tile = tiling[(t_y * this->m_tiling_dims.x) + t_x];
            map[i].sum += tile.sum;
        }

        Float p_x = this->guiding_map.get(i);
        this->m_total_sum += p_x * tile_area;
    }
}

void TileCoding::calc_cdf()
{
    const Point2i dims = this->guiding_map.dims;

    this->m_cdf = std::vector<Float>(dims.y);
    this->m_integral = 0.0;

    float total_pdf_sum = 0.0f;
    for (int y = 0; y < dims.y; ++y)
    {
        float row_pdf_sum = 0.0f;
        for (int x = 0; x < dims.x; ++x)
        {
            int tile_i = (y * dims.x) + x;
            auto p_x = this->guiding_map.get(tile_i) / this->m_total_sum;
            row_pdf_sum += p_x;
        }

        total_pdf_sum += row_pdf_sum;
        this->m_cdf[y] = row_pdf_sum / dims.x;
    }

    this->m_integral = total_pdf_sum / (dims.x * dims.y);
}

void TileCoding::empty_tilings()
{
    this->tilings = std::vector<Tiling>(this->m_tiling_count);
    for (auto& tiling : this->tilings)
    {
        tiling.resize(this->m_tiling_dims.x * this->m_tiling_dims.y);
    }
}

template <typename T>
T TileCoding::transform(T value, bool inverse, const std::function<T(T)>& f)
{
    if (f) return f(value);
    return TCParams::f<T>(this->m_transform_mode, inverse)(value);
}

MTS_NAMESPACE_END