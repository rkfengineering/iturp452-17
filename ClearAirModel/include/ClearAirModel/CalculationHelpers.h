#ifndef INV_CUM_NORM_H
#define INV_CUM_NORM_H

//TODO figure out where to put this. maybe add to mathhelpers

namespace ClearAirModel::CalculationHelpers{

    //From ITU-R P.452-17 Annex 1 Attachment 3

    /// @brief Approximation to inverse cumulative normal distribution function for x<0.5
    /// @param prob threshold probability (fraction)
    /// @return Value from normal distribution (mu = 0, sigma = 1)
    double inv_cum_norm(double prob);

    /// @brief convert frequency (GHz) to wavelength (m) using speed of light constant
    /// @param freq_GHz Frequency (GHz)
    /// @return Wavelength (m)
    double convert_freqGHz_to_wavelength_m(const double& freq_GHz);
}
#endif /* INV_CUM_NORM_H */