#ifndef DIFFRACTION_LOSS_H
#define DIFFRACTION_LOSS_H

#include "ProfilePath.h"
#include "Enumerations.h"

namespace DiffractionLoss {

    /// @brief Bullington part of the diffraction loss from Section 4.2.1
    /// @param d        Vector of distances di of the i-th profile point (km)
    /// @param h        Vector of heights hi of the i-th profile point (masl)
    /// @param hts      Tx Antenna height (masl)
    /// @param hrs      Rx Antenna height (masl)
    /// @param ap       Effective Earth radius (km)
    /// @param freqGHz  Frequency (GHz)
    /// @return Loss from Bullington component (dB)
    double bullLoss(const std::vector<double> d, const std::vector<double> h, const double hts,
            const double hrs, const double ap, const double freqGHz);

    /// @brief Delta-Bullington diffraction loss model from Section 4.2.3
    /// @param path     Contains distance (km) and height (masl) profile points
    /// @param hts      Tx Antenna height (masl)
    /// @param hrs      Rx Antenna height (masl)
    /// @param hstd     Effective height of interfering antenna (masl)
    /// @param hsrd     Effective height of interfered with antenna antenna (masl)
    /// @param ap       Effective Earth radius (km)
    /// @param freqGHz  Frequency (GHz)
    /// @param omega    Fraction of the path over sea
    /// @param pol      Polarization type (horizontal or vertical)
    /// @return Diffraction loss from complete delta-bullington model (dB)
    double delta_bullington(const ProfilePath& path, const double hts, const double hrs, const double hstd, 
            const double hsrd, const double ap, const double freqGHz,
            const double omega, const Enumerations::PolarizationType pol);

    /// @brief Results for the Diffraction Loss model from Section 4.5.4
    /// @param Ldp      Diffraction Loss for the general path not exceeded for p% of the time (dB)
    /// @param Ld50     Diffraction Loss for p=50% (dB)
    struct DiffResults{
        double Ldp;     
        double Ld50;
    };

    /// @brief Diffraction Loss model from Section 4.5.4
    /// @param path     Contains distance (km) and height (masl) profile points
    /// @param hts      Tx Antenna height (masl)
    /// @param hrs      Rx Antenna height (masl)
    /// @param hstd     Effective height of interfering antenna (masl)
    /// @param hsrd     Effective height of interfered with antenna antenna (masl)
    /// @param freqGHz  Frequency (GHz)
    /// @param omega    Fraction of the path over sea
    /// @param pol      Polarization type (horizontal or vertical)
    /// @param p        Percentage of time (%), 0<p<=50
    /// @param b0       Time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100 m of the lower atmosphere (%)
    /// @param DN       Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (positive value) 
    /// @param pol      Polarization type (horizontal or vertical)
    /// @return 
    DiffResults diffLoss(const ProfilePath& path, const double hts, const double hrs, const double hstd, 
            const double hsrd, double freqGHz, const double omega, const double p, const double b0, 
            const double DN, const Enumerations::PolarizationType pol);

    /// @brief Spherical Earth Diffraction Loss exceeded for p% time from Section 4.2.2
    /// @param d_gc     Great circle path distance (km)
    /// @param hte      Effective height of interfering antenna (m)
    /// @param hre      Effective height of interfered-with antenna (m)
    /// @param ap       Effective Earth radius (km)
    /// @param freqGHz  Frequency (GHz)
    /// @param omega    Fraction of path over sea
    /// @param pol      Polarization Type (Horizontal or Vertical)
    /// @return         Spherical-Earth diffraction loss (dB)
    double se_diffLoss(const double d_gc, const double hte, const double hre, const double ap,
            const double freqGHz, const double omega, const Enumerations::PolarizationType pol);


    //maybe hide these in the cpp file (don't expose, keep as local static function)
    
    /// @brief First Term part of Spherical Earth Diffraction Loss from Section 4.2.2.1
    /// @param d_gc     Great circle path distance (km)
    /// @param hte      Effective height of interfering antenna (m)
    /// @param hre      Effective height of interfered-with antenna (m)
    /// @param adft     Effective Earth radius (km)
    /// @param freqGHz  Frequency (GHz)
    /// @param omega    Fraction of path over sea
    /// @param pol      Polarization Type (Horizontal or Vertical)
    /// @return First Term part of Spherical Earth Diffraction Loss (dB)
    double se_first_term(const double d_gc, const double hte, const double hre, const double adft,
            const double freqGHz, const double omega, const Enumerations::PolarizationType pol);
    
    /// @brief Helper Function for First Term part of Spherical Earth Diffraction Loss from Section 4.2.2.1
    /// @param epsr     Relative permittivity 
    /// @param sigma    Conductivity (S/m)
    /// @param d_gc     Great circle path distance (km)
    /// @param hte      Effective height of interfering antenna (m)
    /// @param hre      Effective height of interfered-with antenna (m)
    /// @param adft     Effective Earth radius (km)
    /// @param freqGHz  Frequency (GHz)
    /// @param pol      Polarization Type (Horizontal or Vertical)
    /// @return First Term spherical diffraction loss over a single zone type (dB)
    double se_first_term_inner(const double epsr, const double sigma, const double d_gc, const double hte, const double hre, 
            const double adft, const double freqGHz, const Enumerations::PolarizationType pol);            

}//end namespace DiffractionLoss

#endif /* DIFFRACTION_LOSS_H */