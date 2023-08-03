#ifndef ANOMOLOUS_PROP_H
#define ANOMOLOUS_PROP_H

#include "gtest/gtest.h"
#include "PathProfile.h"
#include "ClearAirModelHelpers.h"

namespace ClearAirModel{

//Section 4.4 Prediction of the basic transmission loss, Lba (dB) 
//occurring during periods of anomalous propagation (ducting and layer reflection)
class AnomalousProp {
    //Allow tests to access private methods
    FRIEND_TEST(MixedProfileTests, AnomalousProp_calcSmoothEarthTxRxHeights_DuctingModel_Test);
    FRIEND_TEST(MixedProfileTests, AnomalousProp_calcTerrainRoughnessTest);

public:
    /// @brief Load inputs for Anomalous Propagation Model and calculate
    /// @param path                     Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param freq_GHz                 Frequency (GHz)
    /// @param height_tx_asl_m          Tx Antenna height (asl_m)
    /// @param height_rx_asl_m          Rx Antenna height (asl_m)
    /// @param temp_K                   Temperature (K)
    /// @param dryPressure_hPa          Dry air pressure (hPa)
    /// @param dist_coast_tx_km         Distance over land from Tx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param dist_coast_rx_km         Distance over land from Rx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param p_percent                Annual percentage of time not exceeded
    /// @param b0_percent               Time percentage that the refractivity gradient exceeds 100 N-Units/km
    /// @param eff_radius_med_km        Median effective Earth's radius (km)
    /// @param horizonVals              Tx and Rx Horizon Elevation Angles (mrad) and Tx and Rx Horizon Distances (km)
    /// @param frac_over_sea            Fraction of the path over sea
    AnomalousProp(const PathProfile::Path& path, const double& freq_GHz,
        const double& height_tx_asl_m, const double& height_rx_asl_m,
        const double& temp_K, const double& dryPressure_hPa, const double& dist_coast_tx_km,
        const double& dist_coast_rx_km, const double& p_percent,
        const double& b0_percent, const double& eff_radius_med_km, 
        const ClearAirModel::HorizonAnglesAndDistances& horizonVals,
        const double& frac_over_sea
    );

    /// @brief get calculated loss value
    /// @return Transmission Loss with ducting and layer reflection (dB)
    inline double getAnomalousPropLoss_dB() const{ return m_anomalousPropLoss_dB; }

private:
    //direct inputs
    const PathProfile::Path& m_path; //Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    const double& m_freq_GHz;        //Frequency (GHz)
    const double& m_height_tx_asl_m; //Tx Antenna height (asl_m)
    const double& m_height_rx_asl_m; //Rx Antenna height (asl_m)
    const double& m_temp_K;          //Temperature (K)
    const double& m_dryPressure_hPa; //Dry air pressure (hPa)
    const double& m_dist_coast_tx_km;//Distance over land from Tx to the coast along the profile path (km) (0 for terminal at sea)
    const double& m_dist_coast_rx_km;//Distance over land from Rx to the coast along the profile path (km) (0 for terminal at sea)
    const double& m_p_percent;       //Annual percentage of time not exceeded

    const double& m_b0_percent;         //Time percentage that the refractivity gradient exceeds 100 N-Units/km
    const double& m_eff_radius_med_km;  //Median effective Earth's radius (km)
    const ClearAirModel::HorizonAnglesAndDistances& m_horizonVals; //Tx and Rx Horizon Elevation Angles (mrad) and Tx and Rx Horizon Distances (km)

    //Consider making this a data member of the path class that gets calculated once
    const double& m_frac_over_sea;      //Fraction of the path over sea

    //calculated internally, exposed for better debugging
    double m_d_tot_km;                              //Distance between Tx and Rx antennas (km)
    ClearAirModel::TxRxPair m_effHeights_ducting_m;//Effective height of Tx and Rx antennas used in ducting/layer reflection model (m)
    double m_terrainRoughness_m;                    //Terrain roughness parameter (m)
    double m_longestContiguousInlandDistance_km;    //Longest contiguous Inland segment in profile path (km)

    //Total Fixed Coupling Losses (except clutter losses) between Antennas and Anomalous propagation structures in atmosphere
    double m_fixedCouplingLoss_dB;  
    //Time percentage and angular-distance dependent losses within the anomalous propagation mechanism
    double m_timePercentageAndAngularDistanceLoss_dB;
    //Basic transmission loss during ducting and layer reflection 
    double m_anomalousPropLoss_dB;

    /// @brief Calculate basic transmission loss during ducting and layer reflection 
    /// @return Transmission Loss (dB)
    double calcAnomalousPropLoss() const;

    /// @brief Total of Fixed Coupling Losses (except clutter losses) between Antennas and Anomalous propagation structures in atmosphere
    /// @return Aggregate Coupling Losses (dB)
    double calcFixedCouplingLoss_helper_dB() const;

    /// @brief Time percentage and angular-distance dependent losses within the anomalous propagation mechanism
    /// @return time percentage and angular-distance dependent losses (dB)
    double calcTimePercentageAndAngularDistanceLoss_helper_dB() const;

    /// @brief Annex 2 Section 5.1.6.4 Calculates effective Antenna Heights for use in the Anomalous Propagation model
    /// @param path Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @return Tx,Rx effective antenna heights for the ducting/layer reflection model (amsl) (m)
    ClearAirModel::TxRxPair calcSmoothEarthTxRxHeights_DuctingModel_amsl_m() const;

    /// @brief Annex 2 Section 5.1.6.4 Calculates maximum height of the terrain above the smooth-Earth 
    ///        surface in the section of the path between, and including, the horizon points
    /// @return Terrain Roughness Parameter (m)
    double calcTerrainRoughness_m() const;

}; //end class AnomalousProp

} //end namespace ClearAirModel
#endif /* ANOMOLOUS_PROP_H */



