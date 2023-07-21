#include "ClearAirModel/AnomolousProp.h"
#include "ClearAirModel/BasicProp.h"
#include "Common/MathHelpers.h"
#include <cmath>
#include <iostream>

double AnomolousProp::calcFixedCouplingLoss_helper_dB(const double& freq_GHz, const EffectiveEarth::HorizonAnglesAndDistances& horizonVals,
        const double& dist_coast_tx_km, const double& dist_coast_rx_km, const double& height_tx_asl_m, 
        const double& height_rx_asl_m, const double& frac_over_sea){

    //Equation 47a Empirical correction to account for the increasing attenuation with wavelength inducted propagation 
    double Alf = 0.0;
    if (freq_GHz<0.5){
        Alf = 45.375 - 137.0*freq_GHz+92.5*freq_GHz*freq_GHz;
    }

    //Equation 48a modified angles for site-shielding diffraction
    const auto [HorizonAngles, HorizonDistances] = horizonVals;
    const auto [horizonElevation_tx_mrad, horizonElevation_rx_mrad] = HorizonAngles;
    const auto [horizonDist_tx_km, horizonDist_rx_km] = HorizonDistances;
    const double mod_horizonElevation_tx_mrad = horizonElevation_tx_mrad - 0.1*horizonDist_tx_km;
    const double mod_horizonElevation_rx_mrad = horizonElevation_rx_mrad - 0.1*horizonDist_rx_km;

    //Equation 48 site-shielding diffraction losses for the interfering station
    double Ast = 0.0;
    if (mod_horizonElevation_tx_mrad>0.0){
        Ast = 20.0*std::log10(1.0+0.361*mod_horizonElevation_tx_mrad*std::sqrt(freq_GHz*horizonDist_tx_km))
            +0.264*mod_horizonElevation_tx_mrad*std::pow(freq_GHz,1.0/3.0);
    }
    //Equation 48 site-shielding diffraction losses for the interfered-with station
    double Asr = 0.0;
    if (mod_horizonElevation_rx_mrad>0.0){
        Asr = 20.0*std::log10(1.0+0.361*mod_horizonElevation_rx_mrad*std::sqrt(freq_GHz*horizonDist_rx_km))
            +0.264*mod_horizonElevation_rx_mrad*std::pow(freq_GHz,1.0/3.0);
    }

    //Equation 49 over-sea surface duct coupling corrections for the interfering station
    double Act = 0.0;
    //Equation 49 over-sea surface duct coupling corrections for the interfered with station
    double Acr = 0.0;
    const bool condition = frac_over_sea>=0.75 
                        && dist_coast_tx_km <= horizonDist_tx_km 
                        && dist_coast_tx_km<= 5.0
                        && dist_coast_rx_km <= horizonDist_rx_km 
                        && dist_coast_rx_km<= 5.0;
    if(condition){
        Act = -3.0 * std::exp(-0.25*dist_coast_tx_km*dist_coast_tx_km)
                * (1.0 + std::tanh(0.07*(50.0-height_tx_asl_m)));
        Acr = -3.0 * std::exp(-0.25*dist_coast_rx_km*dist_coast_rx_km)
                * (1.0 + std::tanh(0.07*(50.0-height_rx_asl_m)));
    }

    //Equation 47
    return 102.45 + 20.0*std::log10(freq_GHz*(horizonDist_tx_km+horizonDist_rx_km))
            + Alf + Ast + Asr + Act + Acr;
}

double AnomolousProp::calcTimePercentageAndAngularDistanceLoss_helper_dB(const double& d_tot_km, const double& freq_GHz,
        const EffectiveEarth::HorizonAnglesAndDistances& horizonVals, const double& eff_radius_med_km, 
        const EffectiveEarth::TxRxPair& effHeights_ducting_m, const double& terrainRoughness_m, 
        const double& longestContiguousInlandDistance_km,  const double& p_percent, const double& b0_percent){

    //Equation 51
    const double specificAttenuation_dB_per_mrad = 5.0e-5*eff_radius_med_km*std::pow(freq_GHz,1.0/3.0);

    //Equation 52a modified angles to remove site shielding component
    const auto [HorizonAngles, HorizonDistances] = horizonVals;
    const auto [horizonElevation_tx_mrad, horizonElevation_rx_mrad] = HorizonAngles;
    const auto [horizonDist_tx_km, horizonDist_rx_km] = HorizonDistances;
    const double corrected_horizonElevation_tx_mrad = std::min(horizonElevation_tx_mrad, 0.1*horizonDist_tx_km);
    const double corrected_horizonElevation_rx_mrad = std::min(horizonElevation_rx_mrad, 0.1*horizonDist_rx_km);

    const double pathAngularDistance_mrad = EffectiveEarth::calcPathAngularDistance_mrad(
        EffectiveEarth::TxRxPair{corrected_horizonElevation_tx_mrad, corrected_horizonElevation_rx_mrad},
        d_tot_km,
        eff_radius_med_km
    );

    //Equation 3a
    const double tau = 1.0 - std::exp(-(4.12e-4*std::pow(longestContiguousInlandDistance_km,2.41)));
    //Equation 55a, epsilon = 3.5
    double alpha = -0.6 - 3.5e-9*std::pow(d_tot_km,3.1)*tau;
    //alpha has a lower limit of -3.4
    alpha = std::max(-3.4, alpha);

    const auto [eff_height_tx_m, eff_height_rx_m] = effHeights_ducting_m;
    //Equation 55 correction for path geometry (mu2)
    const double val1 = MathHelpers::simpleSquare(d_tot_km/(std::sqrt(eff_height_tx_m)+std::sqrt(eff_height_rx_m)));
    const double pathGeometryCorrection = std::pow(500.0/eff_radius_med_km * val1, alpha);

    //Equation 56a Distance beyond horizons of tx and rx, value is limited to at most 40 km
    const double dI = std::min(d_tot_km-horizonDist_tx_km-horizonDist_rx_km, 40.0);
    //Equation 56 correction for terrain roughness (mu3)
    double terrainRoughnessCorrection=1;
    if(terrainRoughness_m>10){
        terrainRoughnessCorrection = std::exp(-4.6e-5*(terrainRoughness_m-10)*(43.0 + 6.0*dI));
    }

    //Equation 54
    const double beta_percent = b0_percent*pathGeometryCorrection*terrainRoughnessCorrection;

    const double log_beta = std::log10(beta_percent);
    //Equation 53a
    const double val2 = -(9.51-4.8*log_beta+0.198*log_beta*log_beta)*1e-6*std::pow(d_tot_km,1.13);
    const double gamma = 1.076/std::pow(2.0058-log_beta,1.012) * std::exp(val2);
    //Equation 53
    const double timePercentageVariabilityLoss_dB = -12.0 + (1.2 + 3.7e-3*d_tot_km)* std::log10(p_percent/beta_percent)
                    + 12.0 * std::pow(p_percent/beta_percent, gamma);
   
    //Equation 50
    return specificAttenuation_dB_per_mrad * pathAngularDistance_mrad + timePercentageVariabilityLoss_dB;
}

double AnomolousProp::calcAnomolousPropLoss(const double& d_tot_km, const double& freq_GHz, const double& height_tx_asl_m,
        const double& height_rx_asl_m, const double& frac_over_sea, const double& temp_K, const double& dryPressure_hPa,
        const double& fixedCouplingLoss_dB, const double& timePercentageAndAngularDistanceLoss_dB){

    //The extra path length from considering the antenna heights is insignificant 
    //but it is also an explicit difference between 452-16 and 452-17
    //Equation 8a distance accounting for height differential (calculated for transhorizon paths too)
    const double d_los_km = std::sqrt(MathHelpers::simpleSquare(d_tot_km)+
                                        MathHelpers::simpleSquare((height_tx_asl_m-height_rx_asl_m)/1000.0));

    //Equation 9a Water Vapor Density
    const double rho = 7.5 + 2.5 * frac_over_sea;
    // Equation 9
    const double gasLoss_dB = BasicProp::calcGasAtten_dB(d_los_km,freq_GHz,temp_K,dryPressure_hPa,rho);

    //Equation 46
    return fixedCouplingLoss_dB+timePercentageAndAngularDistanceLoss_dB+gasLoss_dB;
}
