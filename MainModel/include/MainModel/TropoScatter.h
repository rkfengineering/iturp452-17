#ifndef TROPO_SCATTER_H
#define TROPO_SCATTER_H

#include "PathProfile.h"
#include "ClearAirModelHelpers.h"
#include "DataGridTxt.h"
#include "Common/GeodeticCoord.h"

//TODO figure out if this should be a class

//Section 4.3 Emperical Tropospheric scatter model including secondary propagation effects
// suitable for 0.001% <= p <= 50% (ITU-R P452-17 Section 4.3 NOTE 1)
namespace ITUR_P452::TropoScatter {

    /// @brief Basic transmission loss due to troposcatter not exceeded for p percentage of time
    /// @param commonInputs     //Contains frequency, time percent, antenna heights, and path-related inputs
    /// @param elevationAngles_mrad         Horizon Elevation Angles for transhorizon path, 
    /// @param eff_radius_med_km            Median effective Earth's radius (km)
    /// @param path_angular_distance_mrad   Path Angular Distance (mrad)
    /// @param seaLevelSurfaceRefractivity  Sea Level Surface Refractivity (N0) (N-Units)
    /// @param txHorizonGain_dBi            Tx Antenna directional gain towards the horizon along the path (dB)
    /// @param rxHorizonGain_dBi            Rx Antenna directional gain towards the horizon along the path (dB)
    /// @param temp_K                       Temperature (K)
    /// @param dryPressure_hPa              Dry air pressure (hPa)
    /// @return Loss due to troposcatter (dB)
    double calcTroposcatterLoss_dB(const CommonInputs& commonInputs, const ITUR_P452::TxRxPair&elevationAngles_mrad, 
            const double& eff_radius_med_km, const double& seaLevelSurfaceRefractivity, const double& txHorizonGain_dBi, 
            const double& rxHorizonGain_dBi, const double& temp_K, const double& dryPressure_hPa);
            
}//end namespace TropoScatter
 //end namespace ITUR_P452

#endif /* TROPO_SCATTER_H */