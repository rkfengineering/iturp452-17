#include "BasicProp.h"
#include "GasAttenuationHelpers.h"
#include "MathHelpers.h"
#include <cmath>

double BasicProp::pathLossWithGasAndMultipath(const double& d_tot_km, const double& height_tx_masl, const double& height_rx_masl,
                        const double& freq_GHz, const double& temp_K, const double& dryPressure_hPa, const double& frac_over_sea,
                        const double& d_horizon_t_km, const double& d_horizon_r_km, const double& p_percent){
    
    //Equation 8a distance accounting for height differential
    const double d_los_km = std::sqrt(MathHelpers::simpleSquare(d_tot_km)+
                                        MathHelpers::simpleSquare((height_tx_masl-height_rx_masl)/1000.0));

    //Equation 9a Water Vapor Density
    const double rho = 7.5 + 2.5 * frac_over_sea;
    //Equation 4 from ITU-R P.676-13
    const double waterVapor_hPa = rho * temp_K / 216.7;

    const double totalPressure_hPa = dryPressure_hPa + waterVapor_hPa;
	const double specificAttenuation_dBPerKm = GasAttenuationHelpers::calculateSpecificTotalAttenuation_dBPerKm(
                                            freq_GHz, temp_K, totalPressure_hPa, waterVapor_hPa);
	// Equation 9
	const double gasLoss_dB = specificAttenuation_dBPerKm * d_los_km;

    //Equation 8
    const double freeSpaceWithGasLoss_dB = BasicProp::freeSpacePathLoss(d_los_km,freq_GHz) + gasLoss_dB;

    //Equation 11,12
    return freeSpaceWithGasLoss_dB + BasicProp::multipathFocusingCorrection(d_horizon_t_km, d_horizon_r_km, p_percent);
}

double BasicProp::freeSpacePathLoss(const double& d_los_km, const double& freq_GHz){
    //Equation 8 without gas atten. Modified to avoid calling log10 twice
    return 92.4 + 20.0*std::log10(freq_GHz*d_los_km);
}

double BasicProp::multipathFocusingCorrection(const double& d_horizon_t_km, const double& d_horizon_r_km, const double& p_percent){
    //Eq 10a, 10b 
    return 2.6 * (1.0 - std::exp(-0.1*(d_horizon_t_km+d_horizon_r_km)))*std::log10(p_percent/50.0);
}


