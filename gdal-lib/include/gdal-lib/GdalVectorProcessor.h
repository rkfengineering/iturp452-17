#ifndef GDAL_VECTOR_PROCESSOR_H
#define GDAL_VECTOR_PROCESSOR_H

#include <filesystem>
#include <iostream>
#include <sstream>
#include <string> 
#include <vector>

#include "ogrsf_frmts.h"

class GdalVectorProcessor {
public:
    GdalVectorProcessor(const std::string& vectorFilePath);

    bool doesVectorContainPoint(const double &pixelCenterLat_deg, const double &pixelCenterLon_deg, const double &pixelSize_deg = 0.0) const;

private:
    void processVectorFile();

    std::string vectorFilePath_;
    std::vector<OGRPolygon> vectorPolygonList_;
};

#endif // GDAL_VECTOR_PROCESSOR_H