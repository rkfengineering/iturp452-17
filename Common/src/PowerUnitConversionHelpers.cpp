#include "PowerUnitConversionHelpers.h"

#include <limits>
#include <math.h>

namespace {
    double constexpr min {-std::numeric_limits<double>::max()};
    double constexpr max {std::numeric_limits<double>::max()};
    double constexpr none {-std::numeric_limits<double>::infinity()};
    double constexpr inf {std::numeric_limits<double>::infinity()};
}

const double PowerUnitConversionHelpers::convertDbToWatts(const double& dbValue) {
    // Watts = 10 ^ (dB / 10)
    return std::pow(10.0, dbValue / 10.0);
}

const double PowerUnitConversionHelpers::convertWattsToDb(const double& wattValue) {
    // dB = 10 * log10(Watts)
    return 10.0 * std::log10(wattValue);
}