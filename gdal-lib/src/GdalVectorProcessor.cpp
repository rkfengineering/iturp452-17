#include "gdal-lib/GdalVectorProcessor.h"

GdalVectorProcessor::GdalVectorProcessor(const std::string& vectorFilePath) : vectorFilePath_(vectorFilePath) {
    std::filesystem::path fsPath(vectorFilePath_);
    const bool vectorFileExists = std::filesystem::exists(fsPath);
    
    if (!vectorFileExists) {
        std::ostringstream oStr;
        oStr << "ERROR: GdalVectorProcessor: For reference, the current directory is: " << std::filesystem::current_path() << 
                    "\nThere is no file located at: " << vectorFilePath_ << std::endl;
        throw std::runtime_error(oStr.str());
    }
    processVectorFile();
}

void GdalVectorProcessor::processVectorFile() {
    GDALAllRegister();

    GDALDatasetUniquePtr poDataset(GDALDataset::Open(vectorFilePath_.c_str(), GDAL_OF_VECTOR));
    const bool vectorOpened = poDataset != NULL;

    if (!vectorOpened) {
        std::ostringstream oStr;
        oStr << "ERROR: GdalVectorProcessor: Something went wrong when trying to open the vector file located at: " << vectorFilePath_ << std::endl;
        throw std::runtime_error(oStr.str());
    }

    // Process feature data and store all polygons found in the vector file
    for(OGRLayer* poLayer: poDataset->GetLayers() )
    {
        for( const auto& poFeature: *poLayer )
        {
            const OGRGeometry *poGeometry = poFeature->GetGeometryRef();
            if( poGeometry != nullptr && poGeometry->getGeometryType() == wkbPolygon) {
                const OGRPolygon polygon = *poGeometry->toPolygon();
                vectorPolygonList_.push_back(polygon);
            }
            if( poGeometry != nullptr && poGeometry->getGeometryType() == wkbMultiPolygon) {
                const OGRMultiPolygon multiPolygon = *poGeometry->toMultiPolygon();
                for (auto& polygon : multiPolygon) {
                    vectorPolygonList_.push_back(*polygon);
                }
            }
        }
    }
}

bool GdalVectorProcessor::doesVectorContainPoint(const double &pixelCenterLat_deg, const double &pixelCenterLon_deg, const double &pixelSize_deg) const{
    std::vector<OGRPoint> testPointList;
    
    // Assume that the given test point is in the center of the pixel
    if (pixelSize_deg > 0.0) {
        OGRPoint testPointNorth(pixelCenterLon_deg, pixelCenterLat_deg + pixelSize_deg / 2.0);
        OGRPoint testPointSouth(pixelCenterLon_deg, pixelCenterLat_deg - pixelSize_deg / 2.0);
        OGRPoint testPointEast(pixelCenterLon_deg + pixelSize_deg / 2.0, pixelCenterLat_deg);
        OGRPoint testPointWest(pixelCenterLon_deg - pixelSize_deg / 2.0, pixelCenterLat_deg);
        testPointList = {testPointNorth, testPointSouth, testPointWest, testPointEast};
    }
    // Append the test point itself too
    OGRPoint testPoint(pixelCenterLon_deg, pixelCenterLat_deg);
    testPointList.push_back(testPoint);

    for (const auto& vectorPolygon : vectorPolygonList_) {
        for (const auto& point : testPointList) {
            const bool containsPoint = vectorPolygon.Contains(&point);
            if (containsPoint) {
                return true;
            }
        }
        
    }

    return false;
}