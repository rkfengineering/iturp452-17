#ifndef DIFFRACTION_LOSS_H
#define DIFFRACTION_LOSS_H

#include "gtest/gtest.h"
#include "PathProfile.h"
#include "ClearAirModelHelpers.h"
#include "Common/Enumerations.h"

namespace ITUR_P452{

//Section 4.2 Delta Bullington Diffraction Loss not exceeded for a given annual percentage time
class DiffractionLoss {

FRIEND_TEST(DeltaBullingtonTests, calculateDiffractionModelSmoothEarthHeightsTest);
FRIEND_TEST(DeltaBullingtonTests, DiffractionLoss_calcBullingtonLossTest);
FRIEND_TEST(DeltaBullingtonTests, DiffractionLoss_calcSphericalEarthLossFirstTermTest);
FRIEND_TEST(DeltaBullingtonTests, DiffractionLoss_calcSphericalEarthLossTest);
FRIEND_TEST(DeltaBullingtonTests, DiffractionLoss_calcDeltaBullingtonLossTest);
FRIEND_TEST(DeltaBullingtonTests, DiffractionLoss_calcSmoothEarthBullingtonLossTest);
FRIEND_TEST(MixedProfileTests,DiffractionLossTests_calcSphericalEarthDiffractionLossTest);

public:
    /// @brief Diffraction Loss model from Section 4.5.4
    /// @param path             Contains distance (km) and height (asl)(m) profile points
    /// @param height_tx_asl_m  Tx Antenna height (m)
    /// @param height_rx_asl_m  Rx Antenna height (m)
    /// @param freq_GHz         Frequency (GHz)
    /// @param deltaN           Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (positive value) 
    /// @param pol              Polarization type (horizontal or vertical)
    /// @param p_percent        Percentage of time not exceeded (%), 0<p<=50
    /// @param b0_percent       Time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100 m of the lower atmosphere (%)
    /// @param frac_over_sea    Fraction of the path over sea
    DiffractionLoss(const PathProfile::Path& path, const double& height_tx_asl_m, const double& height_rx_asl_m,
            const double& freq_GHz, const double& deltaN, const ItuModels::Enumerations::PolarizationType& pol, 
            const double& p_percent, const double&b0_percent, const double& frac_over_sea);

    /// @brief Diffraction Loss model from Section 4.5.4
    /// @param out_diff_loss_median_dB Returns diffraction loss not exceeded for 50 percentof time
    /// @param out_diff_loss_p_percent_dB Returns diffraction loss not exceeded for p percent of time
    void calcDiffractionLoss_dB(double& out_diff_loss_median_dB, double& out_diff_loss_p_percent_dB) const;

private:
    //direct inputs
    const PathProfile::Path& m_path; //Contains distance (km) and height (asl)(m) profile points
    const double& m_height_tx_asl_m; //Tx Antenna height (m)
    const double& m_height_rx_asl_m; //Rx Antenna height (m)
    const double& m_freq_GHz;        //Frequency (GHz)
    const double& m_deltaN;          //Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (positive value) 
    const ItuModels::Enumerations::PolarizationType& m_pol; //Polarization type (horizontal or vertical)
    const double& m_p_percent;       // Percentage of time not exceeded (%), 0<p<=50

    //intermediate inputs
    const double& m_b0_percent;      //Time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100 m of the lower atmosphere (%)
    const double& m_frac_over_sea;   //Fraction of the path over sea

    //Calculated intermediate values
    double m_d_tot_km;                //Total great circle path distance from tx to rx (km)
    double m_eff_height_itx_m;        //Effective height of interfering antenna (m)
    double m_eff_height_irx_m;        //Effective height of interfered-with antenna (m)

    ///WARNING When calculating the diffraction parameter, certain square brackets may render as a floor function.
    //They are supposed to be brackets    
    /// @brief Bullington part of the diffraction loss from Section 4.2.1
    /// @param path             Contains distance (km) and height (asl) (m) profile points
    /// @param height_tx_asl_m  Tx Antenna height (asl) (m)
    /// @param height_rx_asl_m  Rx Antenna height (asl) (m)
    /// @param eff_radius_p_km  Effective Earth radius for time percentage (km)
    /// @return Loss from Bullington component (dB)
    double calcBullingtonLoss_dB(const PathProfile::Path& path, const double& height_tx_asl_m,
                                    const double& height_rx_asl_m, const double& eff_radius_p_km) const;

    /// @brief Delta-Bullington diffraction loss model from Section 4.2.3
    /// @param eff_radius_p_km      Effective Earth radius for time percentage (km)
    /// @return Diffraction loss from complete delta-bullington model (dB)
    double calcDeltaBullingtonLoss_dB(const double& eff_radius_p_km) const;

    /// @brief Spherical Earth Diffraction Loss exceeded for p% time from Section 4.2.2
    /// @param eff_radius_p_km       Effective Earth radius for time percentage (km)
    /// @return Spherical-Earth diffraction loss (dB)
    double calcSphericalEarthDiffractionLoss_dB(const double& eff_radius_km) const;
    
    /// @brief First Term part of Spherical Earth Diffraction Loss from Section 4.2.2.1
    /// @param eff_radius_km         Effective Earth radius (km)
    /// @return First Term part of Spherical Earth Diffraction Loss (dB)
    double calcSphericalEarthDiffraction_firstTerm_dB(const double& eff_radius_km) const;
    
    /// @brief Helper Function for First Term part of Spherical Earth Diffraction Loss from Section 4.2.2.1
    /// @param relPermittivity    Relative permittivity 
    /// @param conductivity       Conductivity (S/m)
    /// @param eff_radius_km      Effective Earth radius (km)
    /// @return First Term spherical diffraction loss over a single zone type (dB)
    double calcSphericalEarthDiffraction_firstTerm_singleZone_dB(const double& relPermittivity, const double& conductivity, 
                                                                const double& eff_radius_km) const;

    /// @brief Annex 2 Section 5.1.6.3 Calculates effective Antenna Heights for use in the smooth path Bullington Loss calculation
    ///        in the Delta-Bullington model
    /// @return Tx,Rx effective smooth earth path heights for the diffraction model (amsl) (m)
    ITUR_P452::TxRxPair calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m() const;

};//end class DiffractionLoss
} //end namespace ITUR_P452
#endif /* DIFFRACTION_LOSS_H */


//The Delta Bullington model tries to account for spherical earth diffraction and Bullington diffraction losses in the same model
//If Spherical Earth Diffraction Loss is greater than the bullington loss for the equivalent smooth earth m_path, 
//add the difference to the actual bullington loss. 

//Note: If the m_path is smooth, the bullington losses cancel and the spherical loss term dominates this expression
//For more information on the Delta Bullington model, see 
//https://erdc-library.erdc.dren.mil/jspui/bitstream/11681/42780/1/ERDC-CRREL%20TR-22-1.pdf