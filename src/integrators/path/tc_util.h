#pragma once

#include <cmath>
#include <vector>
#include <random>
#include <array>
#include <functional>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <atomic>
#include <iostream>

#define TC_NAMESPACE_BEGIN namespace TC {
#define TC_NAMESPACE_END }

TC_NAMESPACE_BEGIN

#define TC_Epsilon    1e-08f

#ifndef M_PI
#define M_PI          3.14159265358979323846f
#endif
#ifndef INV_PI
#define INV_PI        0.31830988618379067154f
#endif
#ifndef INV_TWOPI
#define INV_TWOPI     0.15915494309189533577f
#endif
#ifndef INV_FOURPI
#define INV_FOURPI    0.07957747154594766788f
#endif

#define TC_Assert assert

// [[likely]] and [[unlikely]] tag equivalents for C++14
#if defined(__GNUC__) || defined(__clang__)
#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define LIKELY(x)   (x)
#define UNLIKELY(x) (x)
#endif

using Float = float;

struct Sample {
    float phi = 0.0f;
    float theta = 0.0f;
    float value = 0.0f;
    float pdf = 0.0f;
};

template <typename T, typename std::enable_if<std::is_arithmetic<T>::value>::type* = nullptr>
struct Pair {
    T x;
    T y;

    Pair() : x(0), y(0) {};
    Pair(T v) : x(v), y(v) {};
    Pair(T x, T y) : x(x), y(y) {};
};

using Point2 = Pair<float>;
using Point2i = Pair<int>;

struct Util {
private:
    template<typename T>
    class Mapper {
    public:
        explicit Mapper(T value) : value(value) {};

        Mapper<T>& from(const std::pair<T, T>& from_range) {
            this->from_range = from_range;
            return *this;
        }

        T to(const std::pair<T, T>& to_range) {
            if (this->from_range.first == this->from_range.second) {
                return to_range.first;
            }
            T slope = (to_range.second - to_range.first) / (T)(this->from_range.second - this->from_range.first);
            return (slope * (this->value - this->from_range.first) + to_range.first);
        }

    private:
        T value;
        std::pair<T, T> from_range;
    };

public:
    template <typename T, typename std::enable_if<std::is_arithmetic<T>::value>::type* = nullptr>
    static T clamp(const T& n, const T& lower, const T& upper) {
        return std::max(lower, std::min(n, upper));
    };

    static Point2 spherical_to_uv(const Point2& spherical) {
        float u = spherical.x * INV_TWOPI;
        float v = spherical.y * INV_PI;

        return Point2(
            Util::clamp<float>(u, 0.0f, 1.0f - TC_Epsilon),
            Util::clamp<float>(v, 0.0f, 1.0f - TC_Epsilon)
        );
    }

    static Point2 uv_to_spherical(const Point2& uv) {
        float phi = 2.0f * (float) M_PI * uv.x;
        float theta = (float) M_PI * uv.y;
        
        return Point2(
            Util::clamp<float>(phi, 0.0f, 2.0f * (float) M_PI - TC_Epsilon), 
            Util::clamp<float>(theta, 0.0f, (float) M_PI - TC_Epsilon)
        );
    }

    template<typename T>
    static Mapper<T> map(T value) {
        return Mapper<T>(value);
    }
};

struct DomainTransform {
    enum Type {
        Planar,
        Spherical,
        Cosine
    };

    template <typename T>
    static inline T apply(Type t, T value, bool inverse = false) {
        if (t == Planar) {
            return value;
        }
        if (t == Spherical) {
            auto i = std::floor(value);
            auto f = value - i;
            if (inverse) {
                return (1.0 - std::cos(f * M_PI)) * 0.5 + i;
            } else {
                auto pos = (2.0 * f) - 1.0;
                return (std::acos(-pos) * INV_PI) + i;
            }
        }
        if (t == Cosine) {
            auto i = std::floor(value);
            if (inverse) {
                auto f = value - i;
                auto pos = 1.0 - (2.0 * f);
                return (std::acos(pos) * INV_PI) + i;
            } else {
                auto sgn = ((int)i % 2 == 0) ? 1 : -1;
                return 0.5 * (1.0 - sgn * std::cos(value * M_PI)) + i;
            }
        }
        return value;
    }
};

namespace Lookup {
    enum CI : int { P999, P995, P990, P975, P950 };

    struct TTable {
    public:
        static float fetch(CI ci, int n) {
            TC_Assert(n > 0);
            if (n <= 20) return fetch_exact_small(ci, n);

            static const std::array<float, 5> Z_SCORES = { 3.29053f, 2.80703f, 2.57583f, 2.24140f, 1.95996f };
            float z = Z_SCORES[ci];
            float n_f = static_cast<float>(n);
            
            float z2 = z * z; float z3 = z2 * z; float z5 = z3 * z2; float z7 = z5 * z2;
            
            float term1 = (z3 + z) / (4.0f * n_f);
            float term2 = (5.0f * z5 + 16.0f * z3 + 3.0f * z) / (96.0f * n_f * n_f);
            float term3 = (3.0f * z7 + 19.0f * z5 + 17.0f * z3 - 15.0f * z) / (384.0f * n_f * n_f * n_f);
            return z + term1 + term2 + term3;
        }

    private:
        static float fetch_exact_small(CI ci, int n) {
            static const float EXACT[5][20] = {
                { 636.619f, 31.599f, 12.923f, 8.610f, 6.868f, 5.958f, 5.407f, 5.041f, 4.780f, 4.586f, 4.436f, 4.317f, 4.220f, 4.140f, 4.072f, 4.014f, 3.965f, 3.921f, 3.883f, 3.849f },
                { 127.321f, 14.089f, 7.453f, 5.597f, 4.773f, 4.316f, 4.029f, 3.832f, 3.689f, 3.581f, 3.496f, 3.428f, 3.372f, 3.325f, 3.286f, 3.251f, 3.222f, 3.196f, 3.173f, 3.153f },
                { 63.656f, 9.924f, 5.840f, 4.604f, 4.032f, 3.707f, 3.499f, 3.355f, 3.249f, 3.169f, 3.105f, 3.054f, 3.012f, 2.976f, 2.946f, 2.920f, 2.898f, 2.878f, 2.860f, 2.845f },
                { 12.706f, 4.302f, 3.182f, 2.776f, 2.570f, 2.446f, 2.364f, 2.306f, 2.262f, 2.228f, 2.200f, 2.178f, 2.160f, 2.144f, 2.131f, 2.119f, 2.109f, 2.100f, 2.093f, 2.085f },
                { 6.313f, 2.919f, 2.353f, 2.131f, 2.015f, 1.943f, 1.894f, 1.859f, 1.833f, 1.812f, 1.795f, 1.782f, 1.770f, 1.761f, 1.753f, 1.745f, 1.739f, 1.734f, 1.729f, 1.724f }
            };
            return EXACT[ci][n - 1];
        }
    };
}

struct BTCArguments {
    int tilings = 3;
    Point2i tiles = Point2i(1, 1);
    float excess = 0.1f;
    int max_splits = 4095;
    int eagerness = 0; // P999 default
    DomainTransform::Type transformation = DomainTransform::Spherical;
};

// Thread-safe Random Generator
struct RandomGen {
    RandomGen() : m_engine(std::random_device{}()) {}
    RandomGen(uint32_t seed) : m_engine(seed) {}

    float next1D() { return m_distribution(m_engine); }
    Point2 next2D() { return Point2(next1D(), next1D()); }

private:
    std::mt19937 m_engine;
    std::uniform_real_distribution<float> m_distribution{ 0.0f, 1.0f };
};

inline RandomGen& get_thread_random() {
    static thread_local RandomGen generator(std::random_device{}());
    return generator;
}

TC_NAMESPACE_END