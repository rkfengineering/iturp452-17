#ifndef GREAT_CIRCLE_PATH_H
#define GREAT_CIRCLE_PATH_H

#include "gdal-lib/LatLonCoord.h"
#include "gdal-lib/VectorHelpers.h"
#include <vector>
#include <cstdint>

namespace P452 {

class GreatCirclePath{
public:
    GreatCirclePath(const LatLonCoord &startLatLon, const LatLonCoord &endLatLon);

    //interpolates a point a percentage between start and end
    LatLonCoord calcPointAtFractionOfGreatCirclePath_vector(const double& fraction) const;

    std::vector<LatLonCoord> calcPointsOnGreatCirclePath_vector(const uint32_t& numPoints) const;

    /// @brief Calculates points along a great circle path between two given points (start & end points are included in path)
    /// @param startLatLon Starting point of the great circle path
    /// @param endLatLon Ending point of the great circle path
    /// @return A list containing the start, end, and intermediate points along the great circle path between the start & end
    std::vector<LatLonCoord> calcPointsOnGreatCirclePath_sphere(const uint32_t& numPoints) const;

private:
    LatLonCoord startLatLon_, endLatLon_;
    VectorHelpers::ECEFCoordinate startEcef_, endEcef_;
    double startEcefMagnitude_m_, endEcefMagnitude_m_;
    double totalAngle_rad_;
    VectorHelpers::ECEFCoordinate normalStartVector_, normalCrossStartVector_; // For use in intermediate calculations
};

}

#endif // P452