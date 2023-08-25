#ifndef GDAL_RASTER_PROCESSOR_H
#define GDAL_RASTER_PROCESSOR_H

#include <filesystem>
#include <iostream>
#include <sstream>
#include <string> 
#include <vector>

#include "LatLonCoord.h"

class GdalRasterProcessor {
public:
    GdalRasterProcessor(const std::string& rasterFilePath);

    LatLonCoord getNwCorner_latLon() const;
    double getPixelSize_deg() const;
    uint64_t getNumLatitudeRows() const;
    uint64_t getNumLongitudeColumns() const;

    bool fetchValueAtLocation(const LatLonCoord& latLonCoord, double& value) const;
    bool fetchValueAtLocation(const uint64_t& latInd, const uint64_t& lonInd, double& value) const;

private:
    void processRasterFile();

    double nwLatitude_deg_;
    double nwLongitude_deg_;
    double pixelSize_deg_;

    uint64_t numLongitudeColumns_;
    uint64_t numLatitudeRows_;

    std::string rasterFilePath_;
    std::vector<double> rasterData_;
};

#endif // GDAL_RASTER_PROCESSOR_H