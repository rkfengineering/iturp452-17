#ifndef DIFFRACTION_LOSS_H
#define DIFFRACTION_LOSS_H

#include "gtest/gtest.h"
#include "PathProfile.h"
#include "EffectiveEarth.h"
#include "Common/Enumerations.h"

class DiffractionLoss {

FRIEND_TEST(DeltaBullingtonTests, calculateDiffractionModelSmoothEarthHeightsTest);
FRIEND_TEST(DeltaBullingtonTests, DiffractionLoss_calcBullingtonLossTest);
FRIEND_TEST(DeltaBullingtonTests, DiffractionLoss_calcSphericalEarthLossFirstTermTest);
FRIEND_TEST(DeltaBullingtonTests, DiffractionLoss_calcSphericalEarthLossTest);
FRIEND_TEST(DeltaBullingtonTests, DiffractionLoss_calcDeltaBullingtonLossTest);
FRIEND_TEST(DeltaBullingtonTests, DiffractionLoss_calcSmoothEarthBullingtonLossTest);
FRIEND_TEST(MixedProfileTests,DiffractionLossTests_calcSphericalEarthDiffractionLossTest);

public:
    DiffractionLoss(const PathProfile::Path& path, const double& height_tx_asl_m, const double& height_rx_asl_m,
            const double& freq_GHz, const double& deltaN, const Enumerations::PolarizationType& pol, 
            const double& p_percent, const double&b0_percent, const double& frac_over_sea);

    inline double getDiffractionLoss_median_dB() const{ return m_diff_loss_median_dB; }

    inline double getDiffractionLoss_p_percent_dB() const{ return m_diff_loss_p_percent_dB; }

private:
    //direct inputs
    const PathProfile::Path& m_path;
    const double& m_height_tx_asl_m; 
    const double& m_height_rx_asl_m;
    const double& m_freq_GHz;
    const double& m_deltaN;
    const Enumerations::PolarizationType& m_pol;
    const double& m_p_percent;

    //intermediate inputs
    const double& m_b0_percent;
    const double& m_frac_over_sea;

    //Calculated intermediate values
    double m_d_tot_km;
    double m_eff_height_itx_m;        //Effective height of interfering antenna (m)
    double m_eff_height_irx_m;        //Effective height of interfered-with antenna (m)

    //Final Calculated Values
    double m_diff_loss_p_percent_dB;  //Diffraction Loss for the general path not exceeded for p% of the time (dB)
    double m_diff_loss_median_dB;     //Diffraction Loss for p=50% (dB)

    /// @brief Diffraction Loss model from Section 4.5.4
    /// @param path                    Contains distance (km) and height (asl_m) profile points
    /// @param height_tx_asl_m             Tx Antenna height (m)
    /// @param height_rx_asl_m             Rx Antenna height (m)
    /// @param freq_GHz                Frequency (GHz)
    /// @param frac_over_sea           Fraction of the path over sea
    /// @param pol                     Polarization type (horizontal or vertical)
    /// @param p_percent               Percentage of time not exceeded (%), 0<p<=50
    /// @param b0_percent              Time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100 m of the lower atmosphere (%)
    /// @param deltaN                      Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (positive value) 
    /// @param pol                     Polarization type (horizontal or vertical)
    /// @return Returns diffraction loss not exceeded for p percent of time
    double calcDiffractionLoss_p_percent_dB() const;
    double calcDiffractionLoss_median_dB() const;

    ///WARNING When calculating the diffraction parameter, certain square brackets may render as a floor function.
    //They are supposed to be brackets    
    /// @brief Bullington part of the diffraction loss from Section 4.2.1
    /// @param path             Contains distance (km) and height (asl) (m) profile points
    /// @param height_tx_asl_m  Tx Antenna height (asl) (m)
    /// @param height_rx_asl_m  Rx Antenna height (asl) (m)
    /// @param eff_radius_p_km  Effective Earth radius for time percentage (km)
    /// @param freq_GHz         Frequency (GHz)
    /// @return Loss from Bullington component (dB)
    double calcBullingtonLoss_dB(const PathProfile::Path& path, const double& height_tx_asl_m,
                                    const double& height_rx_asl_m, const double& eff_radius_p_km) const;

    /// @brief Delta-Bullington diffraction loss model from Section 4.2.3
    /// @param path                 Contains distance (km) and height (asl_m) profile points
    /// @param height_tx_asl_m          Tx Antenna height (m)
    /// @param height_rx_asl_m          Rx Antenna height (m)
    /// @param eff_radius_p_km      Effective Earth radius for time percentage (km)
    /// @param freq_GHz             Frequency (GHz)
    /// @param frac_over_sea        Fraction of the path over sea
    /// @param pol                  Polarization type (horizontal or vertical)
    /// @return Diffraction loss from complete delta-bullington model (dB)
    double calcDeltaBullingtonLoss_dB(const double& eff_radius_p_km) const;

    //TODO find better names for b0, DN

    /// @brief Spherical Earth Diffraction Loss exceeded for p% time from Section 4.2.2
    /// @param distance_gc_km        Great circle path distance (km)
    /// @param eff_height_itx_m      Effective height of interfering antenna (m)
    /// @param eff_height_irx_m      Effective height of interfered-with antenna (m)
    /// @param eff_radius_p_km       Effective Earth radius for time percentage (km)
    /// @param freq_GHz              Frequency (GHz)
    /// @param frac_over_sea         Fraction of path over sea
    /// @param pol                   Polarization Type (Horizontal or Vertical)
    /// @return Spherical-Earth diffraction loss (dB)
    double calcSphericalEarthDiffractionLoss_dB(const double& eff_radius_km) const;
    
    /// @brief First Term part of Spherical Earth Diffraction Loss from Section 4.2.2.1
    /// @param distance_gc_km        Great circle path distance (km)
    /// @param eff_height_itx_m      Effective height of interfering antenna (m)
    /// @param eff_height_irx_m      Effective height of interfered-with antenna (m)
    /// @param eff_radius_km         Effective Earth radius (km)
    /// @param freq_GHz              Frequency (GHz)
    /// @param frac_over_sea         Fraction of path over sea
    /// @param pol                   Polarization Type (Horizontal or Vertical)
    /// @return First Term part of Spherical Earth Diffraction Loss (dB)
    double calcSphericalEarthDiffraction_firstTerm_dB(const double& eff_radius_km) const;
    
    /// @brief Helper Function for First Term part of Spherical Earth Diffraction Loss from Section 4.2.2.1
    /// @param relPermittivity    Relative permittivity 
    /// @param conductivity       Conductivity (S/m)
    /// @param distance_gc_km     Great circle path distance (km)
    /// @param eff_height_itx_m   Effective height of interfering antenna (m)
    /// @param eff_height_irx_m   Effective height of interfered-with antenna (m)
    /// @param eff_radius_km      Effective Earth radius (km)
    /// @param freq_GHz           Frequency (GHz)
    /// @param pol                Polarization Type (Horizontal or Vertical)
    /// @return First Term spherical diffraction loss over a single zone type (dB)
    double calcSphericalEarthDiffraction_firstTerm_helper_dB(const double& relPermittivity, const double& conductivity, 
                                                                const double& eff_radius_km) const;

    /// @brief Annex 2 Section 5.1.6.3 Calculates effective Antenna Heights for use in the smooth path Bullington Loss calculation
    ///        in the Delta-Bullington model
    /// @param path Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param height_tx_asl_m      Tx antenna height above ground level (m)
    /// @param height_rx_asl_m      Rx antenna height above ground level (m)
    /// @return Tx,Rx effective smooth earth path heights for the diffraction model (amsl) (m)
    EffectiveEarth::TxRxPair calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m() const;

};//end namespace DiffractionLoss

#endif /* DIFFRACTION_LOSS_H */


//The Delta Bullington model tries to account for spherical earth diffraction and Bullington diffraction losses in the same model
//If Spherical Earth Diffraction Loss is greater than the bullington loss for the equivalent smooth earth m_path, 
//add the difference to the actual bullington loss. 

//Note: If the m_path is smooth, the bullington losses cancel and the spherical loss term dominates this expression
//For more information on the Delta Bullington model, see 
//https://erdc-library.erdc.dren.mil/jspui/bitstream/11681/42780/1/ERDC-CRREL%20TR-22-1.pdf