#include "ClearAirModel/P452TotalAttenuation.h"
#include "ClearAirModel/CalculationHelpers.h"

double p452_TotalAttenuation::calcP452TotalAttenuation(const double& freq_GHz, const double& p_percent, const PathProfile::Path path, 
        const double& height_tx_m, const double& height_rx_m, const double& centerLatitude_deg, const double& txHorizonGain_dBi, 
        const double& rxHorizonGain_dBi, const Enumerations::PolarizationType& pol, const double& dist_coast_tx_km, 
        const double& dist_coast_rx_km, const double& deltaN, const double& surfaceRefractivity,
        const double& temp_K, const double& dryPressure_hPa,const double& tx_clutter_height_m, const double& rx_clutter_height_m,
        const double& tx_clutter_dist_km, const double& rx_clutter_dist_km){

    //Path Parameters
    const double effEarthRadius_med_km = EffectiveEarth::calcMedianEffectiveRadius_km(deltaN);
    const double fracOverSea = path.calcFracOverSea();
    const double b0_percent = path.calcTimePercentBeta0(centerLatitude_deg);
    const double longestInlandDistance_km = path.calcLongestContiguousInlandDistance_km();

    //Apply height gain model correction from clutter model
    const auto [mod_path, hg_height_tx_m, hg_height_rx_m, tx_clutterLoss_dB, rx_clutterLoss_dB] = 
            ClutterLoss::calcClutterLoss_dB(freq_GHz,path,height_tx_m,height_rx_m,tx_clutter_height_m,
                                        rx_clutter_height_m,tx_clutter_dist_km,rx_clutter_dist_km);
    const double d_tot_km = mod_path.back().d_km;
    const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
    const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;
    //Fj
    const double slopeInterpolationParameter = 
        p452_TotalAttenuation::calcSlopeInterpolationParameter(mod_path,effEarthRadius_med_km,height_tx_asl_m,height_rx_asl_m);

    //Fk
    const double pathBlendingInterpolationParameter = p452_TotalAttenuation::calcPathBlendingInterpolationParameter(d_tot_km);
    
    const auto HorizonVals = EffectiveEarth::calcHorizonAnglesAndDistances(
        mod_path, height_tx_asl_m, height_rx_asl_m, effEarthRadius_med_km, freq_GHz
    );
    const auto [HorizonAngles, HorizonDistances] = HorizonVals;

    const auto [horizonDist_tx_km, horizonDist_rx_km] = HorizonDistances;

    //TODO rename variables to human-readable


    //TODO can refactor to only calculate free space and gas component once  
    //Free Space Path Loss with Gas attenuation and Multipath/Focusing corrections 
        //Equation 8a distance accounting for height differential (calculated for transhorizon paths too)
    const double d_los_km = std::sqrt(MathHelpers::simpleSquare(d_tot_km)+
                                        MathHelpers::simpleSquare((height_tx_asl_m-height_rx_asl_m)/1000.0));
    //Equation 9a Water Vapor Density
    const double rho = 7.5 + 2.5 * fracOverSea;
    // Equation 9
    const double gasLoss_dB = BasicProp::calcGasAtten_dB(d_los_km,freq_GHz,temp_K,dryPressure_hPa,rho);
    //Equation 8
    const double freeSpaceWithGasLoss_dB = BasicProp::calcFreeSpacePathLoss_dB(d_los_km,freq_GHz) + gasLoss_dB;
    //Equation 10a,10b
    const double calcMultipathFocusingCorrection_p_dB = 
            BasicProp::calcMultipathFocusingCorrection_dB(horizonDist_tx_km, horizonDist_rx_km, p_percent);
    const double calcMultipathFocusingCorrection_b0_dB = 
            BasicProp::calcMultipathFocusingCorrection_dB(horizonDist_tx_km, horizonDist_rx_km, b0_percent);

    const double Lbfsg = freeSpaceWithGasLoss_dB;
    const double Lb0p = freeSpaceWithGasLoss_dB + calcMultipathFocusingCorrection_p_dB;
    const double Lb0b = freeSpaceWithGasLoss_dB + calcMultipathFocusingCorrection_b0_dB;  
    
    const auto DiffractionResults = DiffractionLoss::calcDiffractionLoss_dB(mod_path,hg_height_tx_m, hg_height_rx_m,
            freq_GHz,fracOverSea,p_percent,b0_percent,deltaN,pol);
    
    //Equation 43
    const double Lbd50 = Lbfsg + DiffractionResults.diff_loss_p50_dB;
    //Equation 44
    const double Lbd = Lb0p + DiffractionResults.diff_loss_p_dB;

    //Equation 60
    double Lminb0p = Lb0p + (1-fracOverSea)*DiffractionResults.diff_loss_p_dB;
    if(p_percent>=b0_percent){
        const double Fi = CalculationHelpers::inv_cum_norm(p_percent/100.0)/CalculationHelpers::inv_cum_norm(b0_percent/100.0);
        Lminb0p = Lbd50 + (Lb0b + (1-fracOverSea)*DiffractionResults.diff_loss_p_dB)*Fi;
    }

    const auto eff_heights_ducting = EffectiveEarth::calcSmoothEarthTxRxHeights_DuctingModel_amsl_m(mod_path, hg_height_tx_m,hg_height_rx_m);
    const double terrain_roughness = EffectiveEarth::calcTerrainRoughness_m(mod_path, HorizonDistances);
    
    constexpr double eta = 2.5;
    const double fixed_anomolous = AnomolousProp::calcFixedCouplingLoss_helper_dB(freq_GHz,HorizonVals,dist_coast_tx_km,dist_coast_rx_km,
            height_tx_asl_m, height_rx_asl_m, fracOverSea);
    const double time_anomolous = AnomolousProp::calcTimePercentageAndAngularDistanceLoss_helper_dB(d_tot_km,freq_GHz,HorizonVals,effEarthRadius_med_km,
            eff_heights_ducting, terrain_roughness, longestInlandDistance_km, p_percent, b0_percent);
    const double Lba = AnomolousProp::calcAnomolousPropLoss(d_tot_km, freq_GHz, height_tx_asl_m, height_rx_asl_m, fracOverSea, temp_K, dryPressure_hPa,
            fixed_anomolous, time_anomolous);

    //Equation 61
    const double Lminbap = eta*std::log(std::exp(Lba/eta)+std::exp(Lb0p/eta));
    double Lbda = Lbd;
    if(Lminbap <= Lbd){
        Lbda = Lminbap + (Lbd-Lminbap)*pathBlendingInterpolationParameter;
    }

    //Equation 63
    const double Lbam = Lbda + (Lminb0p-Lbda)* slopeInterpolationParameter;

    const double pathAngularDistance = EffectiveEarth::calcPathAngularDistance_mrad(HorizonAngles,d_tot_km,effEarthRadius_med_km);
    const double Lbs = TropoScatter::calcTroposcatterLoss_dB(d_tot_km,freq_GHz,height_tx_asl_m,height_rx_asl_m,pathAngularDistance,
            surfaceRefractivity, txHorizonGain_dBi, rxHorizonGain_dBi, temp_K, dryPressure_hPa, p_percent);

    const double Lb_pol = -5.0 * std::log10(std::pow(10.0, -0.2*Lbs)+ std::pow(10.0, -0.2*Lbam))+ tx_clutterLoss_dB + rx_clutterLoss_dB;

    //polarization already accounted for
    return Lb_pol;
}

double p452_TotalAttenuation::calcSlopeInterpolationParameter(const PathProfile::Path path, const double effEarthRadius_med_km,
        const double& height_tx_asl_m,const double& height_rx_asl_m){

    const double d_tot = path.back().d_km;
    //Effective Earth Curvature
    const double Ce = 1.0/effEarthRadius_med_km;

    //Code reused from bullington diffraction loss
    //Eq 14 get max slope to profile point from tx 
    double max_slope_tx = std::numeric_limits<double>::lowest();
    double slope_tx;
    PathProfile::ProfilePoint pt_tx;
    for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
        pt_tx = *cit;
        slope_tx = (pt_tx.h_asl_m+500*Ce*pt_tx.d_km*(d_tot-pt_tx.d_km)-height_tx_asl_m)/pt_tx.d_km;
        max_slope_tx = std::max(max_slope_tx,slope_tx); 
    }

    //Eq 15 Slope of line from Tx to Rx assuming LOS
    const double slope_tr_los = (height_rx_asl_m-height_tx_asl_m)/d_tot;
    
    //Adjustable parameters
    constexpr double theta = 0.3;
    constexpr double ksi = 0.8;

    //Equation 58 Calculate interpolation factor Fj
    return 1.0 - 0.5*(1.0 + std::tanh(3.0 * ksi * (max_slope_tx-slope_tr_los)/theta));
}

double p452_TotalAttenuation::calcPathBlendingInterpolationParameter(const double& d_tot_km){
    //fixed parameter for distance range of associated blending
    constexpr double dsw = 20;
    //fixed parameter for blending slope at ends of range
    constexpr double kappa = 0.5;

    //Equation 59 Calculate interpolation factor Fk
    return 1.0 - 0.5 * (1.0 + std::tanh(3.0 * kappa * (d_tot_km-dsw)/dsw));
}
