#include "P452/P452.h"

#include <utility>
#include "MainModel/Dataloader.h"
#include "MainModel/P452TotalAttenuation.h"
#include "GasModel/GasAttenuationHelpers.h"
#include "Common/Enumerations.h"
#include "Common/GeodeticCoord.h"


double P452::calculateP452Loss_dB(const double& txHeight_m, const double& rxHeight_m, 
            const std::vector<double>& elevationList_m, const double& stepDistance_km, 
            const double& midpoint_lat_deg, const double& midpoint_lon_deg,
            const double& freq_GHz, const double& timePercent, const int& polariz,
            const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi,
            const ClutterModel::ClutterType& txClutterType, const ClutterModel::ClutterType& rxClutterType){

    //path creation
    PathProfile::Path p452Path;
    double dist_coast_tx_km,dist_coast_rx_km;
    createP452Path(elevationList_m, stepDistance_km, p452Path, dist_coast_tx_km, dist_coast_rx_km);

    //approximate midpoint height (might not be exactly the midpoint if there are an even number of points)
    //we could check if its even and do an average of the two middle points but I don't think that's worth it
    const double midpointHeight_km = elevationList_m[elevationList_m.size()/2]/1000.0;

    //Get temp_K, dryPressure from standard atmosphere sources
    //assuming summer, mid latitude for Kuwait
    double temp_K, totalPressure_hPa, rho_gm3;
    GasAttenuationHelpers::setAtmosphericTermsForUsMidLatitude(midpointHeight_km, 
                temp_K, totalPressure_hPa, rho_gm3, Enumerations::Season::SummerTime);
    //Equation 4 from ITU-R P.676-13
    const double waterVapor_hPa = rho_gm3 * temp_K / 216.7;
    const double dryPressure_hPa = totalPressure_hPa - waterVapor_hPa;

    //get deltaN, N0 (surfaceRefractivity) from data map
    //convert coordinate to itumodels format
    const GeodeticCoord loc = GeodeticCoord(midpoint_lon_deg, midpoint_lat_deg);
    const double deltaN = ITUR_P452::DataLoader::fetchRadioRefractivityIndexLapseRate(loc);
    const double surfaceRefractivity = ITUR_P452::DataLoader::fetchSeaLevelSurfaceRefractivity(loc);

    //convert polarization convention
    Enumerations::PolarizationType pol;
    if(polariz==0){//itm polarization 0 for horizontal
        pol = Enumerations::PolarizationType::HorizontalPolarized;
    } 
    else{//itm polarization 1 for vertical
        pol = Enumerations::PolarizationType::VerticalPolarized;
    }

    //use ITU-R P.452-17
    const auto p452Model = ITUR_P452::TotalClearAirAttenuation(freq_GHz, timePercent, p452Path, 
            txHeight_m, rxHeight_m, midpoint_lat_deg, txHorizonGain_dBi, 
            rxHorizonGain_dBi, pol, dist_coast_tx_km, dist_coast_rx_km, deltaN, surfaceRefractivity,
            temp_K, dryPressure_hPa, txClutterType, rxClutterType);

    return p452Model.calcTotalClearAirAttenuation();
}

void P452::createP452Path(const std::vector<double>& elevationList_m, const double& stepDistance_km,
        PathProfile::Path& out_path, double& out_dist_coast_tx_km, double& out_dist_coast_rx_km){
   
    //Step 1 convert elevation to path
    //create a new path
    PathProfile::Path newPath;

    //assume starting zone is inland or sea
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

    //Step 2 Fill coastal values
    //go front to back to fill coastal values
    double lastSeaLocation_km = -500.0;//big negative value (or use numeric limits lowest)
    for(auto it = newPath.begin(); it<newPath.end(); ++it){
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
    for(auto it = newPath.rbegin(); it<newPath.rend(); ++it){
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

    //Step 3 find distance to coast
    //These loops have early termination conditions and might not be run at all. Hence why they aren't merged with Step 2
    //distance to coast only matters if its less than 5km. Otherwise we can just put 500km as an arbitrary large value
    //500km is used by P452 validation data for paths that are far from the coast
    out_dist_coast_tx_km=500.0;//initial large value
    out_dist_coast_rx_km=500.0;//initial large value

    if(newPath.front().zone==PathProfile::ZoneType::Sea){
        out_dist_coast_tx_km=0.0;
    }
    else{
        //go front to back
        for(auto it = newPath.begin(); it<newPath.end(); ++it){
            if(it->zone==PathProfile::ZoneType::Sea){
                //end of land area reached (coast reached)
                //take half a step towards the land for border location
                out_dist_coast_rx_km = it->d_km - stepDistance_km/2.0;
                break;
            }
        }
        //if the end of the loop is reached, keep default large value of 500km
    }

    if(newPath.back().zone==PathProfile::ZoneType::Sea){
        out_dist_coast_rx_km=0.0;
    }
    else{
        //go back to front
        for(auto it = newPath.rbegin(); it<newPath.rend(); ++it){
            if(it->zone==PathProfile::ZoneType::Sea){
                //end of land area reached (coast reached)
                //take half a step towards the land for border location
                out_dist_coast_rx_km = newPath.back().d_km-(it->d_km + stepDistance_km/2.0);
                break;
            }
        }
        //if the end of the loop is reached, keep default large value of 500km
    }

    out_path = std::move(newPath);
}
