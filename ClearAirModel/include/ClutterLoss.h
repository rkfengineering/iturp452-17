#ifndef CLUTTER_LOSS_H
#define CLUTTER_LOSS_H

#include "PathProfile.h"

//Section 4.5 Additional clutter losses

namespace ClutterLoss {

    struct ClutterLossResults {
        /// distances (km), heights (masl), and zone types of the profile points in the height gain model
        PathProfile::Path path; 
        /// Tx Antenna center height above ground level (m) in the height gain model
        double hg_height_tx_m;
        /// Rx Antenna center height above ground level (m) in the height gain model
        double hg_height_rx_m;
        /// Additional clutter shielding losses at Tx (dB)
        double tx_clutterLoss_dB;
        /// Additional clutter shielding losses at Rx (dB)
        double rx_clutterLoss_dB;
    };
    
    /// @brief Compute the height-gain correlation from Section 4.5.4. See TABLE 4 for Nominal height/distance inputs
    /// @param freq_GHz             Transmitting Frequency (GHz) 
    /// @param path                 distances (km), heights (masl), and zone types of the profile points
    /// @param height_tx_m          Tx Antenna center height above ground level (m)
    /// @param height_rx_m          Rx Antenna center height above ground level (m)
    /// @param tx_clutter_height_m  Nominal clutter height above ground level at tx (m)
    /// @param rx_clutter_height_m  Nominal clutter height above ground level at rx (m)
    /// @param tx_clutter_dist_km   Distance from nominal clutter point to tx antenna (km)
    /// @param rx_clutter_dist_km   Distance from nominal clutter point to rx antenna (km)
    /// @return         See ClutterLossResults 
    ClutterLossResults clutterLoss_corr(const double& freq_GHz, const PathProfile::Path& path, 
            const double& height_tx_m, const double& height_rx_m, const double& tx_clutter_height_m, const double& rx_clutter_height_m,
            const double& tx_clutter_dist_km, const double& rx_clutter_dist_km);

}//end namespace ClutterLoss

#endif /* CLUTTER_LOSS_H */