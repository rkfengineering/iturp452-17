#include "gdal-lib/GdalRasterProcessor.h"
#include "gdal-lib/LatLonCoord.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <memory>

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()

GdalRasterProcessor::GdalRasterProcessor(const std::string& rasterFilePath) : rasterFilePath_(rasterFilePath) {
    std::filesystem::path fsPath(rasterFilePath_);
    const bool rasterFileExists = std::filesystem::exists(fsPath);
    
    if (!rasterFileExists) {
        std::ostringstream oStr;
        oStr << "ERROR: GdalRasterProcessor: For reference, the current directory is: " << std::filesystem::current_path() << 
                    "\nThere is no file located at: " << rasterFilePath_ << std::endl;
        throw std::runtime_error(oStr.str());
    }
    processRasterFile();
}

void GdalRasterProcessor::processRasterFile() {
    GDALAllRegister();
    GDALDatasetUniquePtr poDataset = GDALDatasetUniquePtr(GDALDataset::Open(rasterFilePath_.c_str(), GA_ReadOnly));
    
    const bool rasterOpened = poDataset != NULL;

    if (!rasterOpened) {
        std::ostringstream oStr;
        oStr << "ERROR: GdalRasterProcessor: Something went wrong when trying to open the raster file located at: " << rasterFilePath_ << std::endl;
        throw std::runtime_error(oStr.str());
    }
    
    double adfGeoTransform[6];

    numLongitudeColumns_ = poDataset->GetRasterXSize();
    numLatitudeRows_ = poDataset->GetRasterYSize();
    std::cout << "Raster size = " << numLongitudeColumns_ << " x " << numLatitudeRows_ << " x " << poDataset->GetRasterCount() << std::endl;
    
    const bool gotGeoTransform = poDataset->GetGeoTransform( adfGeoTransform ) == CE_None;

    if (!gotGeoTransform) {
        std::ostringstream oStr;
        oStr << "ERROR: GdalRasterProcessor: Something went wrong when trying to access the geo-transform from the raster file located at: " << rasterFilePath_ << std::endl;
        throw std::runtime_error(oStr.str());
    }

    nwLongitude_deg_ = adfGeoTransform[0];
    nwLatitude_deg_ = adfGeoTransform[3];

    // Ensure that the pixels are squares
    const double pixelLatSize_deg = adfGeoTransform[5];
    const double pixelLonSize_deg = adfGeoTransform[1];
    assert(abs(abs(pixelLonSize_deg) - abs(pixelLatSize_deg)) < 1.0e-6);
    
    pixelSize_deg_ = abs(pixelLatSize_deg);
    std::cout << "Origin = (" << nwLongitude_deg_ << " deg, " << nwLatitude_deg_ << " deg)" << std::endl;
    std::cout << "Pixel size = (" << pixelLonSize_deg << " deg, " << pixelLatSize_deg << " deg)" << std::endl;
    
    auto poBand = poDataset->GetRasterBand( 1 );
    
    // Reserve memory in the raster data vector
    rasterData_.resize(numLongitudeColumns_ * numLatitudeRows_);
    
    
    // Read all data in the raster file
    int nXOffset = 0;
    int nYOffset = 0;
    int nPixelOffset = 0;
    int nLineOffset = 0;
    auto error = poBand->RasterIO(GF_Read, nXOffset, nYOffset, 
                numLongitudeColumns_, numLatitudeRows_, rasterData_.data(), numLongitudeColumns_, numLatitudeRows_, 
                GDT_Float64, nPixelOffset, nLineOffset);

    if (error != CE_None) {
        std::ostringstream oStr;
        oStr << "ERROR: GdalRasterProcessor: Something went wrong when trying to read data from the raster file located at: " << rasterFilePath_ << std::endl;
        throw std::runtime_error(oStr.str());
    }

    std::cout << "Done processing raster file located at: " << rasterFilePath_ << "!" << std::endl;
}

LatLonCoord GdalRasterProcessor::getNwCorner_latLon() const {
    return LatLonCoord(nwLatitude_deg_, nwLongitude_deg_);
}

double GdalRasterProcessor::getPixelSize_deg() const {
    return pixelSize_deg_;
}

uint64_t GdalRasterProcessor::getNumLatitudeRows() const {
    return numLatitudeRows_; 
}

uint64_t GdalRasterProcessor::getNumLongitudeColumns() const {
    return numLongitudeColumns_; 
}

bool GdalRasterProcessor::fetchValueAtLocation(const LatLonCoord& latLonCoord, double& value) const{
    if (rasterData_.empty()) {
        throw std::runtime_error("ERROR: GdalRasterProcessor: Raster data has been processed yet, cannot fetch anything yet!");
    }

    // Assume NW latitude indicates the N edge of the NW corner's pixel
    const uint64_t latInd = floor((nwLatitude_deg_ - latLonCoord.lat_deg_) / pixelSize_deg_);
    // Assume the NW longitude indicates the W edge of the NW corner's pixel
    const uint64_t lonInd = floor((latLonCoord.lon_deg_ - nwLongitude_deg_) / pixelSize_deg_);

    return fetchValueAtLocation(latInd, lonInd, value);
}

bool GdalRasterProcessor::fetchValueAtLocation(const uint64_t& latInd, const uint64_t& lonInd, double& value) const{
    const uint64_t compositeInd = latInd * numLongitudeColumns_ + lonInd;
    
    if (compositeInd >= rasterData_.size()) {
        return false;
    }

    value = rasterData_[compositeInd];
    return true;
}