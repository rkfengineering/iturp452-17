#ifndef ANOMOLOUS_PROP_H
#define ANOMOLOUS_PROP_H

#include "gtest/gtest.h"
#include "PathProfile.h"
#include "ClearAirModelHelpers.h"

namespace ITUR_P452{

//Section 4.4 Prediction of the basic transmission loss, Lba (dB) 
//occurring during periods of anomalous propagation (ducting and layer reflection)
class AnomalousProp {
    //Allow tests to access private methods
    FRIEND_TEST(MixedProfileTests, AnomalousProp_calcSmoothEarthTxRxHeights_DuctingModel_Test);
    FRIEND_TEST(MixedProfileTests, AnomalousProp_calcTerrainRoughnessTest);

public:
    /// @brief Load inputs for Anomalous Propagation Model and calculate
    /// @param commonInputs             Contains frequency, time percent, antenna heights, and path-related inputs
    /// @param temp_K                   Temperature (K)
    /// @param dryPressure_hPa          Dry air pressure (hPa)
    /// @param dist_coast_tx_km         Distance over land from Tx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param dist_coast_rx_km         Distance over land from Rx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param eff_radius_med_km        Median effective Earth's radius (km)
    /// @param horizonVals              Tx and Rx Horizon Elevation Angles (mrad) and Tx and Rx Horizon Distances (km)
    AnomalousProp(const CommonInputs& commonInputs,
        const double& temp_K, const double& dryPressure_hPa, const double& dist_coast_tx_km,
        const double& dist_coast_rx_km, const double& eff_radius_med_km, 
        const ITUR_P452::HorizonAnglesAndDistances& horizonVals
    );

    /// @brief get calculated loss value
    /// @return Transmission Loss with ducting and layer reflection (dB)
    double calcAnomalousPropLoss_dB() const;

private:
    const CommonInputs& m_commonInputs; //see ClearAirModelHelpers.h

    const double& m_temp_K;          //Temperature (K)
    const double& m_dryPressure_hPa; //Dry air pressure (hPa)

    const double& m_dist_coast_tx_km;//Distance over land from Tx to the coast along the profile path (km) (0 for terminal at sea)
    const double& m_dist_coast_rx_km;//Distance over land from Rx to the coast along the profile path (km) (0 for terminal at sea)

    const double& m_eff_radius_med_km;  //Median effective Earth's radius (km)
    const ITUR_P452::HorizonAnglesAndDistances& m_horizonVals; //Tx and Rx Horizon Elevation Angles (mrad) and Tx and Rx Horizon Distances (km)

    /// @brief Calculate basic transmission loss during ducting and layer reflection due to gaseous attenuation
    /// @return Transmission Loss (dB)
    double calcAnomalousPropGasLoss() const;

    /// @brief Total of Fixed Coupling Losses (except clutter losses) between Antennas and Anomalous propagation structures in atmosphere
    /// @return Aggregate Coupling Losses (dB)
    double calcFixedCouplingLoss_helper_dB() const;

    /// @brief Time percentage and angular-distance dependent losses within the anomalous propagation mechanism
    /// @return time percentage and angular-distance dependent losses (dB)
    double calcTimePercentageAndAngularDistanceLoss_helper_dB() const;

    /// @brief Annex 2 Section 5.1.6.4 Calculates effective Antenna Heights for use in the Anomalous Propagation model
    /// @param path Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @return Tx,Rx effective antenna heights for the ducting/layer reflection model (amsl) (m)
    ITUR_P452::TxRxPair calcSmoothEarthTxRxHeights_DuctingModel_amsl_m() const;

    /// @brief Annex 2 Section 5.1.6.4 Calculates maximum height of the terrain above the smooth-Earth 
    ///        surface in the section of the path between, and including, the horizon points
    /// @return Terrain Roughness Parameter (m)
    double calcTerrainRoughness_m() const;

}; //end class AnomalousProp

} //end namespace ITUR_P452
#endif /* ANOMOLOUS_PROP_H */



