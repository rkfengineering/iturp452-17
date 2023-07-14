#ifndef DIFFRACTION_LOSS_H
#define DIFFRACTION_LOSS_H

#include "PathProfile.h"
#include "Enumerations.h"

namespace DiffractionLoss {

    //WARNING the bullLoss function does not implement floor functions in Eq 14,16,18
    //The ITU published MATLAB implementation ignores these as well
    //The validation data fails when these floor functions are added

    /// @brief Bullington part of the diffraction loss from Section 4.2.1
    /// @param path             Contains distance (km) and height (masl) profile points
    /// @param height_tx_masl   Tx Antenna height (masl)
    /// @param height_rx_masl   Rx Antenna height (masl)
    /// @param eff_radius_p_km  Effective Earth radius for time percentage (km)
    /// @param freqGHz          Frequency (GHz)
    /// @return Loss from Bullington component (dB)
    double bullLoss(const PathProfile::Path& path, const double& height_tx_masl,
            const double& height_rx_masl, const double& eff_radius_p_km, const double& freqGHz);

    /// @brief Delta-Bullington diffraction loss model from Section 4.2.3
    /// @param path                 Contains distance (km) and height (masl) profile points
    /// @param height_tx_masl       Tx Antenna height (masl)
    /// @param height_rx_masl       Rx Antenna height (masl)
    /// @param eff_height_itx_masl  Effective height of interfering antenna (masl)
    /// @param eff_height_irx_masl  Effective height of interfered with antenna (masl)
    /// @param eff_radius_p_km      Effective Earth radius for time percentage (km)
    /// @param freqGHz              Frequency (GHz)
    /// @param frac_over_sea        Fraction of the path over sea
    /// @param pol                  Polarization type (horizontal or vertical)
    /// @return Diffraction loss from complete delta-bullington model (dB)
    double delta_bullington(const PathProfile::Path& path, const double& height_tx_masl, const double& height_rx_masl, 
            const double& eff_height_itx_masl, const double& eff_height_irx_masl, const double& eff_radius_p_km, 
            const double& freqGHz, const double& frac_over_sea, const Enumerations::PolarizationType& pol);

    /// @brief Results for the Diffraction Loss model from Section 4.5.4
    /// @param diff_loss_p_dB      Diffraction Loss for the general path not exceeded for p% of the time (dB)
    /// @param diff_loss_p50_dB    Diffraction Loss for p=50% (dB)
    struct DiffResults{
        double diff_loss_p_dB;     
        double diff_loss_p50_dB;
    };

    //TODO find better names for b0, DN

    /// @brief Diffraction Loss model from Section 4.5.4
    /// @param path                    Contains distance (km) and height (masl) profile points
    /// @param height_tx_masl          Tx Antenna height (masl)
    /// @param height_rx_masl          Rx Antenna height (masl)
    /// @param eff_height_itx_masl     Effective height of interfering antenna (masl)
    /// @param eff_height_irx_masl     Effective height of interfered with antenna antenna (masl)
    /// @param freqGHz                 Frequency (GHz)
    /// @param frac_over_sea           Fraction of the path over sea
    /// @param pol                     Polarization type (horizontal or vertical)
    /// @param p_percent               Percentage of time (%), 0<p<=50
    /// @param b0                      Time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100 m of the lower atmosphere (%)
    /// @param DN                      Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (positive value) 
    /// @param pol                     Polarization type (horizontal or vertical)
    /// @return 
    DiffResults diffLoss(const PathProfile::Path& path, const double& height_tx_masl, const double& height_rx_masl,
            const double& eff_height_itx_masl, const double& eff_height_irx_masl, const double& freqGHz, 
            const double& frac_over_sea, const double& p_percent, const double&b0, 
            const double& DN, const Enumerations::PolarizationType& pol);

    /// @brief Spherical Earth Diffraction Loss exceeded for p% time from Section 4.2.2
    /// @param distance_gc_km        Great circle path distance (km)
    /// @param eff_height_itx_m      Effective height of interfering antenna (m)
    /// @param eff_height_irx_m      Effective height of interfered-with antenna (m)
    /// @param eff_radius_p_km       Effective Earth radius for time percentage (km)
    /// @param freqGHz               Frequency (GHz)
    /// @param frac_over_sea         Fraction of path over sea
    /// @param pol                   Polarization Type (Horizontal or Vertical)
    /// @return Spherical-Earth diffraction loss (dB)
    double se_diffLoss(const double& distance_gc_km, const double& eff_height_itx_m, const double& eff_height_irx_m,
            const double& eff_radius_p_km, const double& freqGHz, const double& frac_over_sea,
            const Enumerations::PolarizationType& pol);


    //maybe hide these in the cpp file (don't expose, keep as local static function)
    
    /// @brief First Term part of Spherical Earth Diffraction Loss from Section 4.2.2.1
    /// @param distance_gc_km        Great circle path distance (km)
    /// @param eff_height_itx_m      Effective height of interfering antenna (m)
    /// @param eff_height_irx_m      Effective height of interfered-with antenna (m)
    /// @param eff_radius_km         Effective Earth radius (km)
    /// @param freqGHz               Frequency (GHz)
    /// @param frac_over_sea         Fraction of path over sea
    /// @param pol                   Polarization Type (Horizontal or Vertical)
    /// @return First Term part of Spherical Earth Diffraction Loss (dB)
    double se_first_term(const double& distance_gc_km, const double& eff_height_itx_m, const double& eff_height_irx_m,
            const double& eff_radius_km, const double& freqGHz, const double& frac_over_sea, 
            const Enumerations::PolarizationType& pol);
    
    /// @brief Helper Function for First Term part of Spherical Earth Diffraction Loss from Section 4.2.2.1
    /// @param eps_r              Relative permittivity 
    /// @param sigma              Conductivity (S/m)
    /// @param distance_gc_km     Great circle path distance (km)
    /// @param eff_height_itx_m   Effective height of interfering antenna (m)
    /// @param eff_height_irx_m   Effective height of interfered-with antenna (m)
    /// @param eff_radius_km      Effective Earth radius (km)
    /// @param freqGHz            Frequency (GHz)
    /// @param pol                Polarization Type (Horizontal or Vertical)
    /// @return First Term spherical diffraction loss over a single zone type (dB)
    double se_first_term_inner(const double& eps_r, const double& sigma, const double& distance_gc_km, 
            const double& eff_height_itx_m, const double& eff_height_irx_m, 
            const double& eff_radius_km, const double& freqGHz, const Enumerations::PolarizationType& pol);            

}//end namespace DiffractionLoss

#endif /* DIFFRACTION_LOSS_H */