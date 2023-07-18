#ifndef BASIC_PROP_H
#define BASIC_PROP_H

#include "PathProfile.h"
//Section 4.1 LOS propagation (including short term effects)
//free space, gas atten, multipath/focusing, 

namespace BasicProp {

    // Note 1. For LOS path, the distance is from the antenna to the Bullington point from the diffraction method for 50% time

    /// @brief LOS transmission loss inlcuding short term effects according to Section 4.1
    /// @param d_tot_km             Distance between Tx and Rx antennas, assuming flat earth (km)
    /// @param height_tx_masl       Tx Antenna height (masl)
    /// @param height_rx_masl       Rx Antenna height (masl)
    /// @param freq_GHz              Frequency (GHz)
    /// @param temp_K               Temperature (K)
    /// @param dryPressure_hPa      Dry air pressure (hPa)
    /// @param frac_over_sea        Fraction of the path over sea
    /// @param d_horizon_t_km       Distance from Tx antenna to its horizon (km), see Note 1
    /// @param d_horizon_r_km       Distance from Rx antenna to its horizon (km), see Note 1
    /// @param p_percent            Percentage of time not exceeded (%), 0<p<=50
    /// @return Transmission Loss (dB)
    double pathLossWithGasAndMultipath(const double& d_tot_km, const double& height_tx_masl, const double& height_rx_masl,
                        const double& freq_GHz, const double& temp_K, const double& dryPressure_hPa, const double& frac_over_sea,
                        const double& d_horizon_t_km, const double& d_horizon_r_km, const double& p_percent);
    
    /// @brief Path loss from Free Space Attenuation
    /// @param d_los_km 
    /// @param freq_GHz 
    /// @return Path Loss (dB)
    double freeSpacePathLoss(const double& d_los_km, const double& freq_GHz);

    /// @brief Corrections for multipath and focusing effects for attenuation not ecveeded for time percentage p
    /// @param d_horizon_t_km       Distance from Tx antenna to its horizon (km), see Note 1
    /// @param d_horizon_r_km       Distance from Rx antenna to its horizon (km), see Note 1
    /// @param p_percent            Percentage of time not exceeded (%), 0<p<=50
    /// @return Attenuation (dB)
    double multipathFocusingCorrection(const double& d_horizon_t_km, const double& d_horizon_r_km, const double& p_percent);

}//end namespace BasicProp

#endif /* BASIC_PROP_H */


