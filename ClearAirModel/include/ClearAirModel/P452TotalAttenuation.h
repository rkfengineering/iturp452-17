#ifndef P452_TOTAL_ATTENUATION_H
#define P452_TOTAL_ATTENUATION_H

#include "PathProfile.h"
#include "EffectiveEarth.h"
#include "BasicProp.h"
#include "DiffractionLoss.h"
#include "TropoScatter.h"
#include "AnomolousProp.h"
#include "ClutterLoss.h"
#include "Common/GeodeticCoord.h"
#include "Common/Enumerations.h"
//Section 4.6

namespace p452_TotalAttenuation {
    /// @brief Calculates Basic transmission loss (dB), not exceeded for the required annual percentage time, p,
    /// @param freq_GHz             Frequency (GHz)
    /// @param p_percent            Percentage of time not exceeded (%), 0<p<=50
    /// @param path                 Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param height_tx_m          Tx Antenna height (m)
    /// @param height_rx_m          Rx Antenna height (m)
    /// @param centerLatitude_deg   The latitude of the path center point (deg) 
    /// @param txHorizonGain_dBi    Tx Antenna directional gain towards the horizon along the path (dB)
    /// @param rxHorizonGain_dBi    Rx Antenna directional gain towards the horizon along the path (dB)
    /// @param pol                  Polarization type (horizontal or vertical)
    /// @param dist_coast_tx_km 
    /// @param dist_coast_rx_km 
    /// @param deltaN 
    /// @param surfaceRefractivity 
    /// @param temp_K 
    /// @param dryPressure_hPa 
    /// @param tx_clutter_height_m 
    /// @param rx_clutter_height_m 
    /// @param tx_clutter_dist_km 
    /// @param rx_clutter_dist_km 
    /// @return 
    double calcP452TotalAttenuation(const double& freq_GHz, const double& p_percent, const PathProfile::Path path, 
            const double& height_tx_m, const double& height_rx_m, const double& centerLatitude_deg, const double& txHorizonGain_dBi, 
            const double& rxHorizonGain_dBi, const Enumerations::PolarizationType& pol, const double& dist_coast_tx_km, 
            const double& dist_coast_rx_km, const double& deltaN, const double& surfaceRefractivity,
            const double& temp_K, const double& dryPressure_hPa,const double& tx_clutter_height_m, const double& rx_clutter_height_m,
            const double& tx_clutter_dist_km, const double& rx_clutter_dist_km);

    double calcSlopeInterpolationParameter(const PathProfile::Path path, const double effEarthRadius_med_km,
            const double& height_tx_asl_m,const double& height_rx_asl_m);
            
    double calcPathBlendingInterpolationParameter(const double& d_tot_km);
}//end namespace p452_TotalAttenuation

#endif /* P452_TOTAL_ATTENUATION_H */