#ifndef CLUTTER_LOSS_H
#define CLUTTER_LOSS_H

#include "ProfilePath.h"

namespace ClutterLoss {

    struct ClutterLossResults {
        /// distances (km), heights (masl), and zone types of the profile points in the height gain model
        ProfilePath path; 
        /// Tx Antenna center height above ground level (m) in the height gain model
        double htgc;
        /// Rx Antenna center height above ground level (m) in the height gain model
        double hrgc;
        /// Additional clutter shielding losses at Tx (dB)
        double Aht;
        /// Additional clutter shielding losses at Rx (dB)
        double Ahr;
    };
    
    /// @brief Compute the height-gain correlation from Section 4.5.4. See TABLE 4 for Nominal height/distance inputs
    /// @param freq_GHz Transmitting Frequency (GHz) 
    /// @param path     distances (km), heights (masl), and zone types of the profile points
    /// @param htg      Tx Antenna center height above ground level (m)
    /// @param hrg      Rx Antenna center height above ground level (m)
    /// @param ha_t     Nominal clutter height above ground level at tx (m)
    /// @param ha_r     Nominal clutter height above ground level at rx (m)
    /// @param dk_t     Distance from nominal clutter point to tx antenna (km)
    /// @param dk_r     Distance from nominal clutter point to rx antenna (km)
    /// @return         See ClutterLossResults (dc,hc,zonec,htgc,hrgc,Aht,Ahr)
    ClutterLossResults closs_corr(const double freq_GHz, const ProfilePath& path, 
            const double htg, const double hrg, const double ha_t, const double ha_r,
            const double dk_t, const double dk_r);

}//end namespace ClutterLoss

#endif /* CLUTTER_LOSS_H */