#ifndef TROPO_SCATTER_H
#define TROPO_SCATTER_H

#include "PathProfile.h"
#include "ClearAirModelHelpers.h"
#include "DataGrid2.h"
#include "Common/GeodeticCoord.h"

//TODO figure out if this should be a class

//Section 4.3 Emperical Tropospheric scatter model including secondary propagation effects
// suitable for 0.001% <= p <= 50% (ITU-R P452-17 Section 4.3 NOTE 1)
namespace ClearAirModel::TropoScatter {

    /// @brief Basic transmission loss due to troposcatter not exceeded for p percentage of time
    /// @param d_tot_km                     Distance between Tx and Rx antennas (km)
    /// @param freq_GHz                     Frequency (GHz)
    /// @param height_tx_asl_m              Tx Antenna height (asl_m)
    /// @param height_rx_asl_m              Rx Antenna height (asl_m)
    /// @param elevationAngles_mrad         Horizon Elevation Angles for transhorizon path, 
    /// @param eff_radius_med_km            Median effective Earth's radius (km)
    /// @param path_angular_distance_mrad   Path Angular Distance (mrad)
    /// @param seaLevelSurfaceRefractivity  Sea Level Surface Refractivity (N0) (N-Units)
    /// @param txHorizonGain_dBi            Tx Antenna directional gain towards the horizon along the path (dB)
    /// @param rxHorizonGain_dBi            Rx Antenna directional gain towards the horizon along the path (dB)
    /// @param temp_K                       Temperature (K)
    /// @param dryPressure_hPa              Dry air pressure (hPa)
    /// @param p_percent                    Percentage of time not exceeded (%), 0.001<=p<=50
    /// @return Loss due to troposcatter (dB)
    double calcTroposcatterLoss_dB(const double& d_tot_km, const double& freq_GHz, const double& height_tx_asl_m,
            const double& height_rx_asl_m, const ClearAirModel::TxRxPair&elevationAngles_mrad, const double& eff_radius_med_km, 
            const double& seaLevelSurfaceRefractivity, const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi, 
            const double& temp_K, const double& dryPressure_hPa, const double& p_percent);

    /// @brief Look up sea level surface refractivity at a location using data map provided with ITU-R P.452 and ITU-R P.1812
    /// @param location Desired location in lon,lat coordinates
    /// @return Refractivity value N0 (N-Units)
    double fetchSeaLevelSurfaceRefractivity(const GeodeticCoord& location);
            
}//end namespace TropoScatter
 //end namespace ClearAirModel

#endif /* TROPO_SCATTER_H */