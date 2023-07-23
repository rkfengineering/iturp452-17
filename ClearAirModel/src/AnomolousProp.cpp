#include "ClearAirModel/AnomolousProp.h"
#include "ClearAirModel/BasicProp.h"
#include "Common/MathHelpers.h"
#include <cmath>
#include <iostream>

AnomolousProp::AnomolousProp(const PathProfile::Path& path, const double& freq_GHz,
    const double& height_tx_asl_m, const double& height_rx_asl_m,
    const double& temp_K, const double& dryPressure_hPa, const double& dist_coast_tx_km,
    const double& dist_coast_rx_km, const double& p_percent,
    const double& b0_percent, const double& eff_radius_med_km, 
    const EffectiveEarth::HorizonAnglesAndDistances& horizonVals,
    const double& frac_over_sea): 
    m_path{path}, m_freq_GHz{freq_GHz}, m_height_tx_asl_m{height_tx_asl_m}, m_height_rx_asl_m{height_rx_asl_m},
    m_temp_K{temp_K}, m_dryPressure_hPa{dryPressure_hPa}, m_dist_coast_tx_km{dist_coast_tx_km}, m_dist_coast_rx_km{dist_coast_rx_km},
    m_p_percent{p_percent}, m_b0_percent{b0_percent}, m_eff_radius_med_km{eff_radius_med_km}, m_horizonVals{horizonVals},
    m_frac_over_sea{frac_over_sea} {

    //Do calculations
    m_d_tot_km = m_path.back().d_km;
    m_effHeights_ducting_m = calcSmoothEarthTxRxHeights_DuctingModel_amsl_m();
    m_terrainRoughness_m = calcTerrainRoughness_m();
    m_longestContiguousInlandDistance_km = m_path.calcLongestContiguousInlandDistance_km();

    //second stage of calculations
    m_fixedCouplingLoss_dB = calcFixedCouplingLoss_helper_dB();
    m_timePercentageAndAngularDistanceLoss_dB = calcTimePercentageAndAngularDistanceLoss_helper_dB();

    //Final stage of calculations
    m_anomolousPropLoss_dB = calcAnomolousPropLoss();
}

double AnomolousProp::calcFixedCouplingLoss_helper_dB()const{

    //Equation 47a Empirical correction to account for the increasing attenuation with wavelength inducted propagation 
    double Alf = 0.0;
    if (m_freq_GHz<0.5){
        Alf = 45.375 - 137.0*m_freq_GHz+92.5*m_freq_GHz*m_freq_GHz;
    }

    //Equation 48a modified angles for site-shielding diffraction
    const auto [HorizonAngles, HorizonDistances] = m_horizonVals;
    const auto [horizonElevation_tx_mrad, horizonElevation_rx_mrad] = HorizonAngles;
    const auto [horizonDist_tx_km, horizonDist_rx_km] = HorizonDistances;
    const double mod_horizonElevation_tx_mrad = horizonElevation_tx_mrad - 0.1*horizonDist_tx_km;
    const double mod_horizonElevation_rx_mrad = horizonElevation_rx_mrad - 0.1*horizonDist_rx_km;

    //Equation 48 site-shielding diffraction losses for the interfering station
    double Ast = 0.0;
    if (mod_horizonElevation_tx_mrad>0.0){
        Ast = 20.0*std::log10(1.0+0.361*mod_horizonElevation_tx_mrad*std::sqrt(m_freq_GHz*horizonDist_tx_km))
            +0.264*mod_horizonElevation_tx_mrad*std::pow(m_freq_GHz,1.0/3.0);
    }
    //Equation 48 site-shielding diffraction losses for the interfered-with station
    double Asr = 0.0;
    if (mod_horizonElevation_rx_mrad>0.0){
        Asr = 20.0*std::log10(1.0+0.361*mod_horizonElevation_rx_mrad*std::sqrt(m_freq_GHz*horizonDist_rx_km))
            +0.264*mod_horizonElevation_rx_mrad*std::pow(m_freq_GHz,1.0/3.0);
    }

    //Equation 49 over-sea surface duct coupling corrections for the interfering station
    double Act = 0.0;
    //Equation 49 over-sea surface duct coupling corrections for the interfered with station
    double Acr = 0.0;
    const bool condition = m_frac_over_sea>=0.75 
                        && m_dist_coast_tx_km <= horizonDist_tx_km 
                        && m_dist_coast_tx_km<= 5.0
                        && m_dist_coast_rx_km <= horizonDist_rx_km 
                        && m_dist_coast_rx_km<= 5.0;
    if(condition){
        Act = -3.0 * std::exp(-0.25*m_dist_coast_tx_km*m_dist_coast_tx_km)
                * (1.0 + std::tanh(0.07*(50.0-m_height_tx_asl_m)));
        Acr = -3.0 * std::exp(-0.25*m_dist_coast_rx_km*m_dist_coast_rx_km)
                * (1.0 + std::tanh(0.07*(50.0-m_height_rx_asl_m)));
    }

    //Equation 47
    return 102.45 + 20.0*std::log10(m_freq_GHz*(horizonDist_tx_km+horizonDist_rx_km))
            + Alf + Ast + Asr + Act + Acr;
}

double AnomolousProp::calcTimePercentageAndAngularDistanceLoss_helper_dB()const{

    //Equation 51
    const double specificAttenuation_dB_per_mrad = 5.0e-5*m_eff_radius_med_km*std::pow(m_freq_GHz,1.0/3.0);

    //Equation 52a modified angles to remove site shielding component
    const auto [HorizonAngles, HorizonDistances] = m_horizonVals;
    const auto [horizonElevation_tx_mrad, horizonElevation_rx_mrad] = HorizonAngles;
    const auto [horizonDist_tx_km, horizonDist_rx_km] = HorizonDistances;
    const double corrected_horizonElevation_tx_mrad = std::min(horizonElevation_tx_mrad, 0.1*horizonDist_tx_km);
    const double corrected_horizonElevation_rx_mrad = std::min(horizonElevation_rx_mrad, 0.1*horizonDist_rx_km);

    const double m_pathAngularDistance_mrad = EffectiveEarth::calcPathAngularDistance_mrad(
        EffectiveEarth::TxRxPair{corrected_horizonElevation_tx_mrad, corrected_horizonElevation_rx_mrad},
        m_d_tot_km,
        m_eff_radius_med_km
    );

    //Equation 3a
    const double tau = 1.0 - std::exp(-(4.12e-4*std::pow(m_longestContiguousInlandDistance_km,2.41)));
    //Equation 55a, epsilon = 3.5
    double alpha = -0.6 - 3.5e-9*std::pow(m_d_tot_km,3.1)*tau;
    //alpha has a lower limit of -3.4
    alpha = std::max(-3.4, alpha);

    const auto [eff_height_tx_m, eff_height_rx_m] = m_effHeights_ducting_m;
    //Equation 55 correction for m_path geometry (mu2)
    const double val1 = MathHelpers::simpleSquare(m_d_tot_km/(std::sqrt(eff_height_tx_m)+std::sqrt(eff_height_rx_m)));
    const double m_pathGeometryCorrection = std::pow(500.0/m_eff_radius_med_km * val1, alpha);

    //Equation 56a Distance beyond horizons of tx and rx, value is limited to at most 40 km
    const double dI = std::min(m_d_tot_km-horizonDist_tx_km-horizonDist_rx_km, 40.0);
    //Equation 56 correction for terrain roughness (mu3)
    double terrainRoughnessCorrection=1;
    if(m_terrainRoughness_m>10){
        terrainRoughnessCorrection = std::exp(-4.6e-5*(m_terrainRoughness_m-10)*(43.0 + 6.0*dI));
    }

    //Equation 54
    const double beta_percent = m_b0_percent*m_pathGeometryCorrection*terrainRoughnessCorrection;

    const double log_beta = std::log10(beta_percent);
    //Equation 53a
    const double val2 = -(9.51-4.8*log_beta+0.198*log_beta*log_beta)*1e-6*std::pow(m_d_tot_km,1.13);
    const double gamma = 1.076/std::pow(2.0058-log_beta,1.012) * std::exp(val2);
    //Equation 53
    const double timePercentageVariabilityLoss_dB = -12.0 + (1.2 + 3.7e-3*m_d_tot_km)* std::log10(m_p_percent/beta_percent)
                    + 12.0 * std::pow(m_p_percent/beta_percent, gamma);
   
    //Equation 50
    return specificAttenuation_dB_per_mrad * m_pathAngularDistance_mrad + timePercentageVariabilityLoss_dB;
}

double AnomolousProp::calcAnomolousPropLoss()const{

    //The extra m_path length from considering the antenna heights is insignificant 
    //but it is also an explicit difference between 452-16 and 452-17
    //Equation 8a distance accounting for height differential (calculated for transhorizon m_paths too)
    const double d_los_km = std::sqrt(MathHelpers::simpleSquare(m_d_tot_km)+
                                        MathHelpers::simpleSquare((m_height_tx_asl_m-m_height_rx_asl_m)/1000.0));

    //Equation 9a Water Vapor Density
    const double rho = 7.5 + 2.5 * m_frac_over_sea;
    // Equation 9
    const double gasLoss_dB = BasicProp::calcGasAtten_dB(d_los_km,m_freq_GHz,m_temp_K,m_dryPressure_hPa,rho);

    //Equation 46
    return m_fixedCouplingLoss_dB+m_timePercentageAndAngularDistanceLoss_dB+gasLoss_dB;
}


EffectiveEarth::TxRxPair AnomolousProp::calcSmoothEarthTxRxHeights_DuctingModel_amsl_m()const{

    //Equations 166a, 166b
    //Tx,Rx heights from a least squares smooth m_path
    auto [height_smooth_tx_amsl_m,height_smooth_rx_amsl_m] = 
            EffectiveEarth::calcLeastSquaresSmoothEarthTxRxHeights_helper_amsl_m(m_path);
    //Equation 168 terminal heights must be above ground level
    height_smooth_tx_amsl_m = std::min(height_smooth_tx_amsl_m, m_path.front().h_asl_m);
    height_smooth_rx_amsl_m = std::min(height_smooth_rx_amsl_m, m_path.back().h_asl_m);

    //Equation 170
    const double eff_height_tx_m = m_height_tx_asl_m - height_smooth_tx_amsl_m;
    const double eff_height_rx_m = m_height_rx_asl_m - height_smooth_rx_amsl_m;

    return EffectiveEarth::TxRxPair(eff_height_tx_m,eff_height_rx_m);
}

//TODO refactor code to make more efficient. We don't have to call the helper function this often
double AnomolousProp::calcTerrainRoughness_m()const{
    //Equations 166a, 166b
    //Tx,Rx heights from a least squares smooth m_path
    auto [height_smooth_tx_amsl_m,height_smooth_rx_amsl_m] = 
            EffectiveEarth::calcLeastSquaresSmoothEarthTxRxHeights_helper_amsl_m(m_path);
    //Equation 168 terminal heights must be above ground level
    height_smooth_tx_amsl_m = std::min(height_smooth_tx_amsl_m, m_path.front().h_asl_m);
    height_smooth_rx_amsl_m = std::min(height_smooth_rx_amsl_m, m_path.back().h_asl_m);

    //smooth earth surface slope
    //assume m_path starts at 0 km 
    const double slope = (height_smooth_rx_amsl_m-height_smooth_tx_amsl_m)/m_path.back().d_km;

    //only evaluate section between horizon points
    const auto [tx_horizon_km, rx_horizon_from_rx] = m_horizonVals.second;//only the distances are needed
    const double rx_horizon_km = m_path.back().d_km -rx_horizon_from_rx;

    //Equation 171 calculate terrain roughness above smooth earth m_path
    double terrainRoughness_m = 0; //the parameter can never be negative 
    double heightAboveSmoothm_path;
    for(auto point : m_path){
        if(point.d_km>=tx_horizon_km && point.d_km<=rx_horizon_km){
            heightAboveSmoothm_path = point.h_asl_m-(height_smooth_tx_amsl_m + slope*point.d_km);
            terrainRoughness_m = std::max(terrainRoughness_m, heightAboveSmoothm_path);
        }
    }
    return terrainRoughness_m;
}
