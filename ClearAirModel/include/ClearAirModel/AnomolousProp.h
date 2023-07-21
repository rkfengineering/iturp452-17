#ifndef ANOMOLOUS_PROP_H
#define ANOMOLOUS_PROP_H

#include "PathProfile.h"
#include "EffectiveEarth.h"

//Section 4.4 Prediction of the basic transmission loss, Lba (dB) 
//occurring during periods of anomalous propagation (ducting and layer reflection)

namespace AnomolousProp {
    /// @brief Calculate basic transmission loss during ducting and layer reflection 
    /// @param d_los_km                 Distance between Tx and Rx antennas (km)
    /// @param freq_GHz                 Frequency (GHz)
    /// @param height_tx_asl_m          Tx Antenna height (asl_m)
    /// @param height_rx_asl_m          Rx Antenna height (asl_m)
    /// @param frac_over_sea            Fraction of the path over sea
    /// @param temp_K                   Temperature (K)
    /// @param dryPressure_hPa          Dry air pressure (hPa)
    /// @param fixedCouplingLoss_dB     Total of Fixed Coupling Losses (except clutter losses) between 
    ///                                 Antennas and Anomolous propagation structures in atmosphere
    /// @param timePercentageAndAngularDistanceLoss_dB Time percentage and angular-distance dependent losses 
    ///                                                within the anomalous propagation mechanism
    /// @return Transmission Loss (dB)
    double calcAnomolousPropLoss(const double& d_tot_km, const double& freq_GHz, const double& height_tx_asl_m,
            const double& height_rx_asl_m, const double& frac_over_sea, const double& temp_K, const double& dryPressure_hPa,
            const double& fixedCouplingLoss_dB, const double& timePercentageAndAngularDistanceLoss_dB);

    /// @brief Total of Fixed Coupling Losses (except clutter losses) between Antennas and Anomolous propagation structures in atmosphere
    /// @param freq_GHz             Frequency (GHz)
    /// @param horizonVals          Tx and Rx Horizon Elevation Angles (mrad) and Tx and Rx Horizon Distances (km)
    /// @param dist_coast_tx_km     Distance over land from Tx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param dist_coast_rx_km     Distance over land from Rx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param height_tx_asl_m      Tx Antenna height above mean sea level (m)
    /// @param height_rx_asl_m      Rx Antenna height above mean sea level (m)
    /// @param frac_over_sea        Fraction of the path over sea
    /// @return Aggregate Coupling Losses (dB)
    double calcFixedCouplingLoss_helper_dB(const double& freq_GHz, const EffectiveEarth::HorizonAnglesAndDistances& horizonVals,
            const double& dist_coast_tx_km, const double& dist_coast_rx_km, const double& height_tx_asl_m, 
            const double& height_rx_asl_m, const double& frac_over_sea);

    /// @brief Time percentage and angular-distance dependent losses within the anomalous propagation mechanism
    /// @param d_tot_km                 Distance between Tx and Rx antennas (km)
    /// @param freq_GHz                 Frequency (GHz)
    /// @param horizonVals              Tx and Rx Horizon Elevation Angles (mrad) and Tx and Rx Horizon Distances (km)
    /// @param eff_radius_med_km        Median effective Earth's radius (km)
    /// @param effHeights_ducting_m     Effective height of Tx and Rx antennas used in ducting/layer reflection model (m)
    /// @param terrainRoughness_m       Terrain roughness parameter (m)
    /// @param longestContiguousInlandDistance_km Longest contiguous Inland segment in profile path (km)
    /// @param b0_percent time percentage that the refractivity gradient exceeds 100 N-Units/km
    /// @return time percentage and angular-distance dependent losses (dB)
    double calcTimePercentageAndAngularDistanceLoss_helper_dB(const double& d_tot_km, const double& freq_GHz,
            const EffectiveEarth::HorizonAnglesAndDistances& horizonVals, const double& eff_radius_med_km, 
            const EffectiveEarth::TxRxPair& effHeights_ducting_m, const double& terrainRoughness_m, 
            const double& longestContiguousInlandDistance_km, const double& p_percent, const double& b0_percent);

}//end namespace AnomolousProp

#endif /* ANOMOLOUS_PROP_H */


/*
    
    //Equation 8a distance accounting for height differential
    const double d_los_km = std::sqrt(MathHelpers::simpleSquare(d_tot_km)+
                                        MathHelpers::simpleSquare((height_tx_asl_m-height_rx_asl_m)/1000.0));
    //Equation 9a Water Vapor Density
    const double rho = 7.5 + 2.5 * frac_over_sea;
   
    // Equation 9
    const double gasLoss_dB = BasicProp::gasAttenWrapper_dB(d_los_km,freq_GHz,temp_K,dryPressure_hPa,rho);
    
    */
