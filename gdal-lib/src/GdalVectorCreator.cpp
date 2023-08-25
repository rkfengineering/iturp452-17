#include "gdal-lib/GdalVectorCreator.h"

#include <filesystem>
#include <iostream>
#include <sstream>

#include "ogrsf_frmts.h"

void GdalVectorCreator::createPointCluster_GeoJsonFile(const std::string& geoJsonFilePath, const std::map<LatLonCoord, double> pointToValueMap) {
    const std::string driverName = "GeoJSON";

    GDALAllRegister();

    const auto poDriver = GetGDALDriverManager()->GetDriverByName(driverName.c_str());
    if( poDriver == NULL )
    {
        std::cout << "ERROR: GdalVectorProcessor::createGeoJsonFile: " << driverName << " driver not available." << std::endl;
        return;
    }

    const auto poDataset = poDriver->Create(geoJsonFilePath.c_str(), 0, 0, 0, GDT_Unknown, NULL);
    if( poDataset == NULL )
    {
        std::cout << "ERROR: GdalVectorProcessor::createGeoJsonFile: Failed to create the output file at " << geoJsonFilePath << std::endl;
        return;
    }

    const auto poLayer = poDataset->CreateLayer("PointClusterLayer", NULL, wkbPoint, NULL );
    if( poLayer == NULL )
    {
        std::cout << "ERROR: GdalVectorProcessor::createGeoJsonFile: Failed to create a layer" << std::endl;
        return;
    }

    OGRFieldDefn oField("value", OFTReal);

    if( poLayer->CreateField( &oField ) != OGRERR_NONE )
    {
        std::cout << "ERROR: GdalVectorProcessor::createGeoJsonFile: Failed to create a field" << std::endl;
        return;
    }

    for (const auto &[latLonCoord, value] : pointToValueMap) {
        auto poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
        poFeature->SetField("value", value);

        // NOTE: x = longitude, y = latitude
        OGRPoint latLonPoint;
        latLonPoint.setX(latLonCoord.lon_deg_);
        latLonPoint.setY(latLonCoord.lat_deg_);

        poFeature->SetGeometry(&latLonPoint);

        if(poLayer->CreateFeature(poFeature) != OGRERR_NONE)
        {
            std::cout << "ERROR: GdalVectorProcessor::createGeoJsonFile: Failed to create one of the given points as a feature in the GeoJSON file" << std::endl;
            return;
        }

        OGRFeature::DestroyFeature(poFeature);
    }

    GDALClose(poDataset);
}

std::string GdalVectorCreator::createPolygonWkt(const std::vector<LatLonCoord> openCoordList) {
    std::ostringstream oStrStream;
    oStrStream << "POLYGON ((";
    for (const auto& coord : openCoordList) {
        oStrStream << coord.lon_deg_ << " " << coord.lat_deg_ << ",";
    }
    oStrStream << openCoordList.front().lon_deg_ << " " << openCoordList.front().lat_deg_ << "))";

    return oStrStream.str();
}

void GdalVectorCreator::createPixels_GeoJsonFile(const std::string& geoJsonFilePath, const std::map<LatLonCoord, double> pixelCenterToValueMap, 
            const double& pixelSize_deg) {
    const std::string driverName = "GeoJSON";

    GDALAllRegister();

    const auto poDriver = GetGDALDriverManager()->GetDriverByName(driverName.c_str());
    if( poDriver == NULL )
    {
        std::cout << "ERROR: GdalVectorProcessor::createGeoJsonFile: " << driverName << " driver not available." << std::endl;
        return;
    }

    const auto poDataset = poDriver->Create(geoJsonFilePath.c_str(), 0, 0, 0, GDT_Unknown, NULL);
    if( poDataset == NULL )
    {
        std::cout << "ERROR: GdalVectorProcessor::createGeoJsonFile: Failed to create the output file at " << geoJsonFilePath << std::endl;
        return;
    }

    const auto poLayer = poDataset->CreateLayer("PixelLayer", NULL, wkbPoint, NULL );
    if( poLayer == NULL )
    {
        std::cout << "ERROR: GdalVectorProcessor::createGeoJsonFile: Failed to create a layer" << std::endl;
        return;
    }

    OGRFieldDefn oField("value", OFTReal);

    if( poLayer->CreateField( &oField ) != OGRERR_NONE )
    {
        std::cout << "ERROR: GdalVectorProcessor::createGeoJsonFile: Failed to create a field" << std::endl;
        return;
    }

    const double halfPixelSize_deg = pixelSize_deg / 2.0;
    for (const auto &[pixelCenterCoord, value] : pixelCenterToValueMap) {
        auto poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
        poFeature->SetField("value", value);
        
        const std::vector<LatLonCoord> pixelPolygonCoordList = {
            LatLonCoord(pixelCenterCoord.lat_deg_ + halfPixelSize_deg, pixelCenterCoord.lon_deg_ + halfPixelSize_deg),
            LatLonCoord(pixelCenterCoord.lat_deg_ + halfPixelSize_deg, pixelCenterCoord.lon_deg_ - halfPixelSize_deg),
            LatLonCoord(pixelCenterCoord.lat_deg_ - halfPixelSize_deg, pixelCenterCoord.lon_deg_ - halfPixelSize_deg),
            LatLonCoord(pixelCenterCoord.lat_deg_ - halfPixelSize_deg, pixelCenterCoord.lon_deg_ + halfPixelSize_deg),
        };
        
        std::string pixelWkt = GdalVectorCreator::createPolygonWkt(pixelPolygonCoordList);
        auto pixelWkt_cStr = pixelWkt.c_str();
        OGRPolygon pixelPolygon;
        pixelPolygon.importFromWkt(&pixelWkt_cStr);

        poFeature->SetGeometry(&pixelPolygon);

        if(poLayer->CreateFeature(poFeature) != OGRERR_NONE)
        {
            std::cout << "ERROR: GdalVectorProcessor::createGeoJsonFile: Failed to create one of the given pixel polygons as a feature in the GeoJSON file" << std::endl;
            return;
        }

        OGRFeature::DestroyFeature(poFeature);
    }

    GDALClose(poDataset);
}
