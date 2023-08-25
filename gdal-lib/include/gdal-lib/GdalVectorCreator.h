#ifndef GDAL_VECTOR_CREATOR_H
#define GDAL_VECTOR_CREATOR_H

#include <map>
#include <string>
#include <vector>

#include "LatLonCoord.h"

namespace GdalVectorCreator {
    std::string createPolygonWkt(const std::vector<LatLonCoord> openCoordList);

    void createPointCluster_GeoJsonFile(const std::string& geoJsonFilePath, const std::map<LatLonCoord, double> pointToValueMap);
    void createPixels_GeoJsonFile(const std::string& geoJsonFilePath, const std::map<LatLonCoord, double> pixelCenterToValueMap, 
            const double& pixelSize_deg);
}

#endif // GDAL_VECTOR_CREATOR_H