#ifndef EFFECTIVE_EARTH_H
#define EFFECTIVE_EARTH_H

#include "PathProfile.h"

namespace EffectiveEarth{

    //Pair for returning tx (first) and rx (second) values
    using TxRxPair = std::pair<double,double>;

    //Pair for returning horizon angles(mrad) (first) and horizon distances(km) (second)
    using HorizonAnglesAndDistances = std::pair<TxRxPair,TxRxPair>;

    //Eq 6b, Median Effective Earth radius exceeded for b% of time
    //kb = 3 at Point Incidence of anomalous propagation b0
    constexpr double k_eff_radius_bpercentExceeded_km = 6371.0 * 3.0; 

    /// @brief Median value of effective earth radius exceeded for 50% of time
    /// @param delta_N  Average radio-refractivity lapse-rate through the lowest 1 km of the atmosphere (N-units/km)
    /// @return Effective Earth Radius (km)
    double calcMedianEffectiveRadius_km(const double& delta_N);

    /// @brief Annex 2 Section 5.1.6.2 Creates a smooth least-squares straight line approximation for the path
    ///        Helper method for calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m
    ///        WARNING Does not account for Eq 168
    /// @param path Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @return Tx,Rx endpoint heights for the smooth-earth surface (amsl) (m)
    TxRxPair calcLeastSquaresSmoothEarthTxRxHeights_helper_amsl_m(const PathProfile::Path& path);

    /// @brief Annex 2 Section 5.1.6.3 Calculates effective Antenna Heights for use in the smooth path Bullington Loss calculation
    ///        in the Delta-Bullington model
    /// @param path Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param height_tx_m      Tx antenna height above ground level (m)
    /// @param height_rx_m      Rx antenna height above ground level (m)
    /// @return Tx,Rx effective antenna heights for the diffraction model (amsl) (m)
    TxRxPair calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m(const PathProfile::Path& path, const double& height_tx_m, const double& height_rx_m);

    /// @brief Annex 1 Attachment 2 Section 4,5 Calculating Antenna Horizon Elevation Angle and Horizon Distances
    /// @param path                 Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param height_tx_asl_m      Tx Antenna height (asl_m)
    /// @param height_rx_asl_m      Rx Antenna height (asl_m)
    /// @param eff_radius_med_km    Median effective Earth's radius (km)
    /// @param freq_GHz             Frequency (GHz)
    /// @return Antenna Horizon Distances (km) and Horizon Elevation Angles (mrad)
    HorizonAnglesAndDistances calcHorizonAnglesAndDistances(const PathProfile::Path& path, const double& height_tx_asl_m,
                                const double& height_rx_asl_m, const double& eff_radius_med_km, const double& freq_GHz);
    
    /// @brief Calculates the path angular distance from the path profile analysis results
    /// @param elevationAngles_mrad Horizon Elevation Angles for transhorizon path, 
    ///                             otherwise elevation angles towards other antenna for LOS path
    /// @param dtot_km              Distance between Tx and Rx (km)
    /// @param eff_radius_med_km    Median effective Earth's radius (km)
    /// @return Path Angular Distance parameter (mrad)
    double calcPathAngularDistance_mrad(const TxRxPair& elevationAngles_mrad, const double& dtot_km, const double& eff_radius_med_km);

    /// @brief Annex 2 Section 5.1.6.4 Calculates effective Antenna Heights for use in the Anomolous Propagation model
    /// @param path Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param height_tx_m      Tx antenna height above ground level (m)
    /// @param height_rx_m      Rx antenna height above ground level (m)
    /// @return Tx,Rx effective antenna heights for the ducting/layer reflection model (amsl) (m)
    TxRxPair calcSmoothEarthTxRxHeights_DuctingModel_amsl_m(const PathProfile::Path& path, const double& height_tx_m, const double& height_rx_m);

    /// @brief Annex 2 Section 5.1.6.4 Calculates maximum height of the terrain above the smooth-Earth 
    ///        surface in the section of the path between, and including, the horizon points
    /// @param path Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param horizonDists_km Horizon distance points from Tx and Rx (km)
    /// @return Terrain Roughness Parameter (m)
    double calcTerrainRoughness_m(const PathProfile::Path& path, TxRxPair horizonDists_km);
}
#endif /* EFFECTIVE_EARTH_H */