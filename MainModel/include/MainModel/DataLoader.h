#ifndef ITUR_P452_DATA_LOADER_H
#define ITUR_P452_DATA_LOADER_H

#include "DataGridTxt.h"
#include "Common/GeodeticCoord.h"

#include <memory>

//TODO consider making this a singleton
namespace ITUR_P452 {
    class DataLoader {
    public:
        /// @brief Fetch sea level surface refractivity 
        /// @param location Desired location in lon,lat coordinates
        /// @return Refractivity value N0 (N-Units)
        /// This data is fetched from a TXT file from ITU-R P.452 and ITU-R P.1812
        static double fetchSeaLevelSurfaceRefractivity(const ItuModels::GeodeticCoord& location);

        /// @brief Fetch Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere
        /// @param location Desired location in lon,lat coordinates
        /// @return Refractivity Lapse Rate Delta N (N-Units/km)
        /// This data is fetched from a TXT file from ITU-R P.452 and ITU-R P.1812
        static double fetchRadioRefractivityIndexLapseRate(const ItuModels::GeodeticCoord& location);

    private:
        //sea level surface refractivity (N-Units)
        static const DataGridTxt _seaLevelSurfaceRefractivityMap;
        //Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (N-Units/km)
        static const DataGridTxt _radioRefractivityIndexLapseRateMap;
    };
} // end namespace Gas

#endif /* ITUR_P452_DATA_LOADER_H */