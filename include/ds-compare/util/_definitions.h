// This file contains definitions (e.g., includes, macros, etc.) that are to be included by each util file in this folder.

#pragma once

#if !defined(__DSCOMPARE_UTIL_DEFINITIONS_H_)
#define __DSCOMPARE_UTIL_DEFINITIONS_H_

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/random.h>

#include <Eigen/Core>

#include <boost/algorithm/clamp.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include <random>
#include <array>
#include <stack>
#include <functional>
#include <unordered_map>
#include <ctime>
#include <ratio>
#include <chrono>
#include <bit>
#include <numeric>

#if !defined(DS_COMPARE)
    #define DS_COMPARE MTS_IMPORT
#endif

typedef std::unordered_map<uint32_t, std::vector<mitsuba::Float>> umap_samples;

#endif /* __DSCOMPARE_UTIL_DEFINITIONS_H_ */