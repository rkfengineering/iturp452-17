#include "P452/P452.h"
#include "P452/GreatCirclePath.h"

#include <utility>
#include "MainModel/Dataloader.h"
#include "MainModel/P452TotalAttenuation.h"
#include "GasModel/GasAttenuationHelpers.h"
#include "Common/Enumerations.h"
#include "Common/GeodeticCoord.h"

//use raw elevation data inputs
double P452::calculateP452Loss_dB(const double& txHeight_m, const double& rxHeight_m, 
            const std::vector<double>& elevationList_m, const double& stepDistance_km, 
            const double& midpoint_lat_deg, const double& midpoint_lon_deg,
            const double& freq_GHz, const double& timePercent, const int& polariz,
            const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi,
            const ClutterModel::ClutterType& txClutterType, const ClutterModel::ClutterType& rxClutterType){

    //P452 path creation from raw elevation
    PathProfile::Path p452Path;
    createP452Path(elevationList_m, stepDistance_km, p452Path);

    //use P452 path
    return calculateP452Loss_dB(txHeight_m, rxHeight_m, p452Path, midpoint_lat_deg, midpoint_lon_deg,
        freq_GHz, timePercent, polariz, txHorizonGain_dBi, rxHorizonGain_dBi, txClutterType, rxClutterType);

}

//use gdal input for elevation data
double P452::calculateP452Loss_dB(const double& txHeight_m, const double& rxHeight_m, 
            const LatLonCoord startCoord, const LatLonCoord endCoord,
            const std::vector<GdalRasterProcessor>& rasterProcessorList,
            const double& freq_GHz, const double& timePercent, const int& polariz,
            const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi,
            const ClutterModel::ClutterType& txClutterType,
            const ClutterModel::ClutterType& rxClutterType){

    //Step 1 create great circle path
    const uint32_t numPoints = 50;
    const GreatCirclePath path(startCoord, endCoord);
    //50 points for now, we can probably use the pixel size and path distance to figure out an appropriate value
    const auto pathPointList = path.calcPointsOnGreatCirclePath_vector(numPoints);
    const LatLonCoord midPointCoord = path.calcPointAtFractionOfGreatCirclePath_vector(0.5);

    //Step 2 create raw elevation list at great circle path points
    const auto startEcef = VectorHelpers::ConvertLatLonToECEF(pathPointList.front());
    const auto nextEcef = VectorHelpers::ConvertLatLonToECEF(pathPointList[1]);
    const double pathStepDistance_km = VectorHelpers::CalculateDistanceMeters(startEcef, nextEcef)/1000.0;

    std::vector<double> rawElevationList;

    // Fill in rest of path elevation list with elevation data
    double elevation_m = 0.0;
    //uint32_t numPointsInSea = 0;
    for (const auto& pathPoint : pathPointList) {
        bool foundElevation = false;
        // Loop through each terrain file until we find the one that contains the location
        for (const auto& rasterProcessor : rasterProcessorList) {
            if (rasterProcessor.fetchValueAtLocation(pathPoint, elevation_m)) {
                rawElevationList.push_back(elevation_m);
                // We fetched the elevation succesfully, go to the next point
                foundElevation = true;
                break;
            }
        }
        // Make sure we still set a path elevation (use 0 as the default)
        if (!foundElevation) {
            rawElevationList.push_back(0.0);
        }
    }

    //P452 path creation from raw elevation
    PathProfile::Path p452Path;
    createP452Path(rawElevationList, pathStepDistance_km, p452Path);

    //use P452 path
    return calculateP452Loss_dB(txHeight_m, rxHeight_m, p452Path, midPointCoord.lat_deg_, midPointCoord.lon_deg_,
        freq_GHz, timePercent, polariz, txHorizonGain_dBi, rxHorizonGain_dBi, txClutterType, rxClutterType);
}


double P452::calculateP452Loss_dB(const double& txHeight_m, const double& rxHeight_m, 
            const LatLonCoord startCoord, const LatLonCoord endCoord,
            const std::vector<GdalRasterProcessor>& rasterProcessorList,
            const GdalVectorProcessor& landBorders,
            const double& freq_GHz, const double& timePercent, const int& polariz,
            const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi,
            const ClutterModel::ClutterType& txClutterType,
            const ClutterModel::ClutterType& rxClutterType){

    PathProfile::Path p452Path;
    LatLonCoord midpointCoord;
    //P452 path creation from raw elevation from gdal inputs
    createP452Path(startCoord, endCoord, rasterProcessorList, landBorders,
                p452Path, midpointCoord);

    //use P452 path
    return calculateP452Loss_dB(txHeight_m, rxHeight_m, p452Path, midpointCoord.lat_deg_, midpointCoord.lon_deg_,
        freq_GHz, timePercent, polariz, txHorizonGain_dBi, rxHorizonGain_dBi, txClutterType, rxClutterType);
}

//wrapper for main entry point in algorithm. Fetches environment values and makes certain assumptions
double P452::calculateP452Loss_dB(const double& txHeight_m, const double& rxHeight_m, 
            const PathProfile::Path& p452path, const double& midpoint_lat_deg, const double& midpoint_lon_deg,
            const double& freq_GHz, const double& timePercent, const int& polariz,
            const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi,
            const ClutterModel::ClutterType& txClutterType, const ClutterModel::ClutterType& rxClutterType){

    //calculate distance to coast
    double dist_coast_tx_km,dist_coast_rx_km;
    calcCoastDistance_km(p452path, dist_coast_tx_km, dist_coast_rx_km);

    //convert coordinate to itumodels format
    double midpointHeight_km;
    uint32_t midIndex = p452path.size()/2;
    if (p452path.size()%2==0){
        midpointHeight_km = (p452path[midIndex].h_asl_m+p452path[midIndex-1].h_asl_m)/2000.0;
    }
    else{
        midpointHeight_km = p452path[midIndex].h_asl_m/1000.0;
    }
    const GeodeticCoord midpointCoord = GeodeticCoord(midpoint_lon_deg, midpoint_lat_deg, midpointHeight_km);

    //get deltaN, N0 (surfaceRefractivity) from data map
    const double deltaN = ITUR_P452::DataLoader::fetchRadioRefractivityIndexLapseRate(midpointCoord);
    const double surfaceRefractivity = ITUR_P452::DataLoader::fetchSeaLevelSurfaceRefractivity(midpointCoord);

    //Get temp_K, dryPressure from standard atmosphere sources
    //assuming summer, mid latitude for Kuwait
    double temp_K, totalPressure_hPa, waterVapor_hPa;

    GasAttenuationHelpers::setSeasonalAtmosphericTermsForUsLocation(midpointCoord, 
                temp_K, totalPressure_hPa, waterVapor_hPa, Enumerations::Season::SummerTime);

    const double dryPressure_hPa = totalPressure_hPa - waterVapor_hPa;

    //convert polarization convention
    Enumerations::PolarizationType pol;
    if(polariz==0){//itm polarization 0 for horizontal
        pol = Enumerations::PolarizationType::HorizontalPolarized;
    } 
    else{//itm polarization 1 for vertical
        pol = Enumerations::PolarizationType::VerticalPolarized;
    }

    //use ITU-R P.452-17
    const auto p452Model = ITUR_P452::TotalClearAirAttenuation(freq_GHz, timePercent, p452path, 
            txHeight_m, rxHeight_m, midpoint_lat_deg, txHorizonGain_dBi, 
            rxHorizonGain_dBi, pol, dist_coast_tx_km, dist_coast_rx_km, deltaN, surfaceRefractivity,
            temp_K, dryPressure_hPa, txClutterType, rxClutterType);

    return p452Model.calcTotalClearAirAttenuation();
}

//create path from raw elevation data
void P452::createP452Path(const std::vector<double>& elevationList_m, const double& stepDistance_km,
        PathProfile::Path& out_path){
   
    //Step 1 convert elevation to path
    //create a new path
    PathProfile::Path newPath;

    //assign points as land or sea
    double distance_km = 0;
    PathProfile::ZoneType zone;
    for (const double& elevation_m : elevationList_m) {
        if(elevation_m==0){
            zone=PathProfile::ZoneType::Sea; //assume sea if elevation is 0
        }
        else{
            zone=PathProfile::ZoneType::Inland;
        }
        newPath.push_back(PathProfile::ProfilePoint(distance_km, elevation_m, zone));
        distance_km+=stepDistance_km;   
    }

    //Step 2 assign coastal zone to qualified inland points
    modifyPathAddCoastalValues(newPath);

    out_path = std::move(newPath);
}

//create path from gdal inputs
void P452::createP452Path(const LatLonCoord startCoord, const LatLonCoord endCoord,
        const std::vector<GdalRasterProcessor>& rasterProcessorList, const GdalVectorProcessor& landBorders,
        PathProfile::Path& out_path, LatLonCoord& out_midpointCoord){
   
    //Step 1 create great circle path
    //create a new path
    PathProfile::Path newPath;

    //50 points for now, we can probably use the pixel size and path distance to figure out an appropriate value
    const uint32_t numPoints = 50;
    const GreatCirclePath path(startCoord, endCoord);
    const auto pathPointList = path.calcPointsOnGreatCirclePath_vector(numPoints);
    out_midpointCoord = path.calcPointAtFractionOfGreatCirclePath_vector(0.5);

    //Step 2 get elevation and zone at great circle path points
    const auto startEcef = VectorHelpers::ConvertLatLonToECEF(pathPointList.front());
    const auto nextEcef = VectorHelpers::ConvertLatLonToECEF(pathPointList[1]);
    const double pathStepDistance_km = VectorHelpers::CalculateDistanceMeters(startEcef, nextEcef)/1000.0;

    std::vector<double> rawElevationList;

    //assume starting zone is inland or sea
    double distance_km = 0.0;
    PathProfile::ZoneType zone;
    double elevation_m = 0.0;
    for (const auto& pathPoint : pathPointList) {
        bool foundElevation = false;
        //use land borders to assign zone
        if(landBorders.doesVectorContainPoint(pathPoint.lat_deg_, pathPoint.lon_deg_)){
            zone=PathProfile::ZoneType::Inland;
        }
        else{
            zone=PathProfile::ZoneType::Sea;
        }

        // Loop through each terrain file until we find the one that contains the location
        for (const auto& rasterProcessor : rasterProcessorList) {
            if (rasterProcessor.fetchValueAtLocation(pathPoint, elevation_m)) {
                newPath.push_back(PathProfile::ProfilePoint(distance_km, elevation_m, zone));
                // We fetched the elevation succesfully, go to the next point
                foundElevation = true;
                break;
            }
        }
        // Make sure we still set a path elevation (use 0 as the default)
        if (!foundElevation) {
            newPath.push_back(PathProfile::ProfilePoint(distance_km, 0.0, zone)); 
        }
        distance_km+=pathStepDistance_km;
    }

    //step 3 add coastal zone data
    modifyPathAddCoastalValues(newPath);

    out_path = std::move(newPath);
}

void P452::calcCoastDistance_km(const PathProfile::Path& path, double& out_dist_coast_tx_km, double& out_dist_coast_rx_km){
    //distance to coast only matters if its less than 5km. Otherwise we can just put 500km as an arbitrary large value
    //500km is used by P452 validation data for paths that are far from the coast
    out_dist_coast_tx_km=500.0;//initial large value
    out_dist_coast_rx_km=500.0;//initial large value

    //assume constant path spacing
    const double stepDistance_km = path.at(1).d_km;

    if(path.front().zone==PathProfile::ZoneType::Sea){
        out_dist_coast_tx_km=0.0;
    }
    else{
        //go front to back
        for(auto it = path.cbegin(); it<path.cend(); ++it){
            if(it->zone==PathProfile::ZoneType::Sea){
                //end of land area reached (coast reached)
                //take half a step towards the land for border location
                out_dist_coast_rx_km = it->d_km - stepDistance_km/2.0;
                break;
            }
        }
        //if the end of the loop is reached, keep default large value of 500km
    }

    if(path.back().zone==PathProfile::ZoneType::Sea){
        out_dist_coast_rx_km=0.0;
    }
    else{
        //go back to front
        for(auto it = path.crbegin(); it<path.crend(); ++it){
            if(it->zone==PathProfile::ZoneType::Sea){
                //end of land area reached (coast reached)
                //take half a step towards the land for border location
                out_dist_coast_rx_km = path.back().d_km-(it->d_km + stepDistance_km/2.0);
                break;
            }
        }
        //if the end of the loop is reached, keep default large value of 500km
    }
}

void P452::modifyPathAddCoastalValues(PathProfile::Path& path){
    //go front to back to fill coastal values
    double lastSeaLocation_km = -500.0;//big negative value (or use numeric limits lowest)
    for(auto it = path.begin(); it<path.end(); ++it){
        if(it->zone==PathProfile::ZoneType::Sea){
            lastSeaLocation_km = it->d_km;
        }
        else if (it->zone==PathProfile::ZoneType::Inland && it->h_asl_m<=100.0){
            //only consider points within 50km of sea
            if(it->d_km-lastSeaLocation_km<=50.0){
                it->zone = PathProfile::ZoneType::CoastalLand;
            }
        }
    }
    //go back to front to fill coastal values
    lastSeaLocation_km = 500.0;//big positive value (or use numeric limits max)
    for(auto it = path.rbegin(); it<path.rend(); ++it){
        if(it->zone==PathProfile::ZoneType::Sea){
            lastSeaLocation_km = it->d_km;
        }
        else if (it->zone==PathProfile::ZoneType::Inland && it->h_asl_m<=100.0){
            //only consider points within 50km of sea
            if(lastSeaLocation_km-it->d_km<=50.0){
                it->zone = PathProfile::ZoneType::CoastalLand;
            }
        }
    }
}
