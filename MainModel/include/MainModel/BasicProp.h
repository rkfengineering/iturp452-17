#ifndef BASIC_PROP_H
#define BASIC_PROP_H

#include "gtest/gtest.h"
#include "PathProfile.h"
#include "ClearAirModelHelpers.h"

namespace ITUR_P452{

//Section 4.1 LOS propagation (including short term effects)
//free space, gas atten, multipath/focusing, 

class BasicProp {
    FRIEND_TEST(BasicPropTests, calcFreeSpacePathLoss);
public:
    /// @brief Load inputs for Basic Propagation Model with gaseous attenuation and multipath focusing and calculate
    /// @param commonInputs     //Contains frequency, time percent, antenna heights, and path-related inputs
    /// @param temp_K           //Temperature (K)
    /// @param dryPressure_hPa  //Dry air pressure (hPa)
    /// @param horizonDists_km  //Tx and Rx Horizon Distances (km)
    BasicProp(const CommonInputs& commonInputs, const double& temp_K, const double& dryPressure_hPa, 
                const ITUR_P452::TxRxPair& horizonDists_km);

    /// @brief Calculate values from the basic transmission loss model of Section 4.1
    /// @param out_freeSpaceWithGasLoss_dB              return free space transmission loss with gas attenuation
    /// @param out_basicTransmissionLoss_p_percent_dB   return free space loss with gas atten and multipath focusing correction for p percent of time
    /// @param out_basicTransmissionLoss_b0_percent_dB  return free space loss with gas atten and multipath focusing correction for b0 percent of time
    void calcTransmissionlosses_dB(double& out_freeSpaceWithGasLoss_dB, double& out_basicTransmissionLoss_p_percent_dB, 
                                    double& out_basicTransmissionLoss_b0_percent_dB) const;

private:
    //direct inputs
    const CommonInputs& m_commonInputs; //see ClearAirModelHelpers.h
    const double& m_temp_K;          //Temperature (K)
    const double& m_dryPressure_hPa; //Dry air pressure (hPa)
    // Note: For LOS path, these distances are from the antenna to the Bullington point from the diffraction method for 50% time
    const ITUR_P452::TxRxPair& m_horizonDists_km; //Tx and Rx Horizon Distances (km)

    /// @brief LOS transmission loss including Gaseous attenuation
    /// @return Transmission Loss (dB)
    double calcPathLossWithGas_dB() const;
    
    /// @brief Path loss from Free Space Attenuation
    /// @param d_los_km Distance (km)
    /// @param freq_GHz Frequency (GHz)
    /// @return Path Loss (dB)
    static double calcFreeSpacePathLoss_dB(const double& d_los_km, const double& freq_GHz);

    /// @brief Corrections for multipath and focusing effects for attenuation not exceeded for time percentage p
    /// @param p_percent            Percentage of time not exceeded (%), 0<p<=50
    /// @return Attenuation (dB)
    double calcMultipathFocusingCorrection_dB(const double& p_percent) const;

};//end class BasicProp
} //end namespace ITUR_P452

#endif /* BASIC_PROP_H */


