#include "ClearAirModel/BasicProp.h"
#include "GasModel/GasAttenuationHelpers.h"
#include "Common/MathHelpers.h"
#include <cmath>

double BasicProp::calcPathLossWithGasAndMultipath_dB(const double& d_tot_km, const double& height_tx_asl_m, const double& height_rx_asl_m,
                        const double& freq_GHz, const double& temp_K, const double& dryPressure_hPa, const double& frac_over_sea,
                        const double& d_horizon_t_km, const double& d_horizon_r_km, const double& p_percent){
    
    //Equation 8a distance accounting for height differential (calculated for transhorizon paths too)
    const double d_los_km = std::sqrt(MathHelpers::simpleSquare(d_tot_km)+
                                        MathHelpers::simpleSquare((height_tx_asl_m-height_rx_asl_m)/1000.0));
    //Equation 9a Water Vapor Density
    const double rho = 7.5 + 2.5 * frac_over_sea;
   
    // Equation 9
    const double gasLoss_dB = BasicProp::calcGasAtten_dB(d_los_km,freq_GHz,temp_K,dryPressure_hPa,rho);

    //Equation 8
    const double freeSpaceWithGasLoss_dB = BasicProp::calcFreeSpacePathLoss_dB(d_los_km,freq_GHz) + gasLoss_dB;

    //Equation 10a,10b
    const double calcMultipathFocusingCorrection_dB = BasicProp::calcMultipathFocusingCorrection_dB(d_horizon_t_km, d_horizon_r_km, p_percent);

    //Equation 11,12
    return freeSpaceWithGasLoss_dB + calcMultipathFocusingCorrection_dB;
}

double BasicProp::calcFreeSpacePathLoss_dB(const double& d_los_km, const double& freq_GHz){
    //Equation 8 without gas atten. Modified to avoid calling log10 twice
    return 92.4 + 20.0*std::log10(freq_GHz*d_los_km);
}

double BasicProp::calcMultipathFocusingCorrection_dB(const double& d_horizon_t_km, const double& d_horizon_r_km, const double& p_percent){
    //Eq 10a, 10b 
    return 2.6 * (1.0 - std::exp(-0.1*(d_horizon_t_km+d_horizon_r_km)))*std::log10(p_percent/50.0);
}

double BasicProp::calcGasAtten_dB(const double& d_los_km, const double& freq_GHz, const double& temp_K, 
                                        const double& dryPressure_hPa, const double& waterVaporDensity_g_m3){

    //Equation 4 from ITU-R P.676-13
    const double waterVapor_hPa = waterVaporDensity_g_m3 * temp_K / 216.7;

    const double totalPressure_hPa = dryPressure_hPa + waterVapor_hPa;

    //Validation data from ITU uses P676-13 with frequencies below 1 GHz 
    const double specificAttenuation_dBPerKm = GasAttenuationHelpers::calculateSpecificTotalAttenuation_dBPerKm(
                                            freq_GHz, temp_K, totalPressure_hPa, waterVapor_hPa);
    // Equation 9
    return specificAttenuation_dBPerKm * d_los_km;
}