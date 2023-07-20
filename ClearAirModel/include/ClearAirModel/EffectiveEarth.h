#ifndef EFFECTIVE_EARTH_H
#define EFFECTIVE_EARTH_H

#include "PathProfile.h"

namespace EffectiveEarth{

    //Eq 6b, Median Effective Earth radius exceeded for b% of time
    //kb = 3 at Point Incidence of anomalous propagation b0
    constexpr double k_eff_radius_bpercentExceeded_km = 6371.0 * 3.0; 

    /// @brief Median value of effective earth radius exceeded for 50% of time
    /// @param delta_N  Average radio-refractivity lapse-rate through the lowest 1 km of the atmosphere (N-units/km)
    /// @return Effective Earth Radius (km)
    double calcMedianEffectiveRadius_km(const double& delta_N);

    //struct for returning tx/rx height pairs
    struct TxRxPair{
        double tx_val;
        double rx_val;
    };

    /// @brief Annex 2 Section 5.1.6.2 Creates a smooth least-squares straight line approximation for the path
    ///        Helper method for calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m
    ///        WARNING Does not account for Eq 168
    /// @param path Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @return Tx,Rx endpoint heights for the smooth-earth surface (amsl) (m)
    TxRxPair calcLeastSquaresSmoothEarthTxRxHeights_helper_amsl_m(const PathProfile::Path& path);

    /// @brief Annex 2 Section 5.1.6.3 Calculates effective Surface Heights for use in the smooth path Bullington Loss calculation
    ///        in the Delta-Bullington model
    /// @param path Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param height_tx_m      Tx antenna height above ground level (m)
    /// @param height_rx_m      Rx antenna height above ground level (m)
    /// @return Tx,Rx effective antenna heights for the diffraction model (amsl) (m)
    TxRxPair calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m(const PathProfile::Path& path, const double& height_tx_m, const double& height_rx_m);

    //Struct for returning horizon angles(mrad) and horizon distances(km)
    //These values are intrinsically connected and would be difficult to calculate separately
    struct HorizonAnglesAndDistances{
        TxRxPair horizonDist_km;
        TxRxPair horizonElevation_mrad;
    };

    /// @brief Annex 1 Attachment 2 Section 4,5 Calculating Antenna Horizon Elevation Angle and Horizon Distances
    /// @param path                 Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param height_tx_asl_m       Tx Antenna height (asl_m)
    /// @param height_rx_asl_m       Rx Antenna height (asl_m)
    /// @param eff_radius_med_km    Median effective Earth's radius (km)
    /// @param freq_GHz             Frequency (GHz)
    /// @return Antenna Horizon Distances (km) and Horizon Elevation Angles (mrad)
    HorizonAnglesAndDistances calcHorizonAnglesAndDistances(const PathProfile::Path& path, const double& height_tx_asl_m,
                                const double& height_rx_asl_m, const double& eff_radius_med_km, const double& freq_GHz);


//TODO implement function for equation 159

}
#endif /* EFFECTIVE_EARTH_H */