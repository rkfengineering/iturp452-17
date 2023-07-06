
#ifndef DECIBEL_H
#define DECIBEL_H

// Helper functions to convert between a decibel value (dBW, dBm) to a linear value (Watts, milliWatts)
namespace PowerUnitConversionHelpers {
    const double convertDbToWatts(const double& dbValue);
    const double convertWattsToDb(const double& wattValue);
};

#endif /* POWER_UNIT_CONVERSION_HELPERS_H */