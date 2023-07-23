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
    
    //Path geometry parameters of modified path
    const auto HorizonVals = EffectiveEarth::calcHorizonAnglesAndDistances(
        mod_path, height_tx_asl_m, height_rx_asl_m, effEarthRadius_med_km, freq_GHz
    );
    const auto [HorizonAngles, HorizonDistances] = HorizonVals;
    const auto [horizonDist_tx_km, horizonDist_rx_km] = HorizonDistances;

    //TODO move this into troposcatter section
    const double pathAngularDistance = EffectiveEarth::calcPathAngularDistance_mrad(HorizonAngles,d_tot_km,effEarthRadius_med_km);

    //Equation 8 (Lbfsg)
    const double freeSpaceWithGasLoss_dB = BasicProp::calcPathLossWithGas_dB(d_tot_km,height_tx_asl_m,
                                            height_rx_asl_m,freq_GHz,temp_K,dryPressure_hPa,fracOverSea);
    //Equation 10a
    const double calcMultipathFocusingCorrection_p_dB = 
            BasicProp::calcMultipathFocusingCorrection_dB(horizonDist_tx_km, horizonDist_rx_km, p_percent);
    //Equation 10b
    const double calcMultipathFocusingCorrection_b0_dB = 
            BasicProp::calcMultipathFocusingCorrection_dB(horizonDist_tx_km, horizonDist_rx_km, b0_percent);

    //TODO rename variables to human-readable

    //Equation 11 (Lb0p) basic transmission loss with gas and multipath not exceeded for p percent of time
    const double basicTransmissionLoss_p_percent_dB = freeSpaceWithGasLoss_dB + calcMultipathFocusingCorrection_p_dB;
    //Equation 12 (Lb0b)
    const double basicTransmissionLoss_b0_percent_dB = freeSpaceWithGasLoss_dB + calcMultipathFocusingCorrection_b0_dB;  
    
    //Delta Bullington Diffraction Loss calculations
    const auto DiffractionModel = DiffractionLoss(
        mod_path, height_tx_asl_m, height_rx_asl_m, freq_GHz, deltaN, pol, p_percent, b0_percent, fracOverSea
    );
    const double diffractionLoss_p_percent_dB = DiffractionModel.getDiffractionLoss_p_percent_dB();
    
    //Equation 43 (Lbd50) basic loss with diffraction loss not exceeded for 50 percent of time
    const double basicWithMedianDiffractionLoss_dB = freeSpaceWithGasLoss_dB + DiffractionModel.getDiffractionLoss_median_dB();
    //Equation 44 (Lbd) basic loss with diffraction loss not exceeded for p percent of time
    const double basicWithDiffractionLoss_p_percent_dB = basicTransmissionLoss_p_percent_dB + diffractionLoss_p_percent_dB;

    //Equation 60 Minimum basic transmission loss associated with LOS propagation and over-sea sub-path diffraction
    double minLossWithOverSeaSubPathDiffraction_dB = basicTransmissionLoss_p_percent_dB + (1-fracOverSea)*diffractionLoss_p_percent_dB;
    if(p_percent>=b0_percent){
        //Diffraction Interpolation parameter
        const double diffractionInterpolationParameter = 
                    CalculationHelpers::inv_cum_norm(p_percent/100.0)/CalculationHelpers::inv_cum_norm(b0_percent/100.0);
                    
        minLossWithOverSeaSubPathDiffraction_dB = MathHelpers::interpolate1D(
                basicWithMedianDiffractionLoss_dB, 
                basicTransmissionLoss_b0_percent_dB + (1-fracOverSea)*diffractionLoss_p_percent_dB,
                diffractionInterpolationParameter
        );
    }

    //Anomolous Propagation Calculations (Ducting and Layer Reflection)
    //TODO figure out a better way to handle inputs/where to put what calculations
    const auto AnomolousPropModel = AnomolousProp(mod_path, freq_GHz, height_tx_asl_m, height_rx_asl_m,
        temp_K, dryPressure_hPa, dist_coast_tx_km, dist_coast_rx_km, p_percent,
        b0_percent, effEarthRadius_med_km, HorizonVals, fracOverSea);

    const double AnomolousPropagationLoss_dB = AnomolousPropModel.getAnomolousPropLoss_dB();

    //Equation 61 (Lminbap)
    constexpr double eta = 2.5; //constant parameter
    const double minLossWithAnomolousPropagation_dB = 
                eta*std::log(std::exp(AnomolousPropagationLoss_dB/eta)+std::exp(basicTransmissionLoss_p_percent_dB/eta));

    //Equation 62 (Lbda)
    double diffractionAndAnomolousPropagationLoss_dB = basicWithDiffractionLoss_p_percent_dB;
    if(minLossWithAnomolousPropagation_dB <= basicWithDiffractionLoss_p_percent_dB){
        diffractionAndAnomolousPropagationLoss_dB = MathHelpers::interpolate1D(minLossWithAnomolousPropagation_dB,
                                                    basicWithDiffractionLoss_p_percent_dB,pathBlendingInterpolationParameter);
    }

    //Equation 63 (Lbam)
    const double modifiedDiffractionAndAnomolousPropagationLoss_dB = 
                                                    MathHelpers::interpolate1D(diffractionAndAnomolousPropagationLoss_dB, 
                                                    minLossWithOverSeaSubPathDiffraction_dB,slopeInterpolationParameter);

    //Calculate Tropospheric Scatter
    const double tropoScatterLoss_dB = TropoScatter::calcTroposcatterLoss_dB(d_tot_km,freq_GHz,height_tx_asl_m,height_rx_asl_m,
            pathAngularDistance, surfaceRefractivity, txHorizonGain_dBi, rxHorizonGain_dBi, temp_K, dryPressure_hPa, p_percent);

    //Equation 64 Total Loss predicted by model, combines losses using a geometric mean of the linear values
    const double val1 = std::pow(10.0, -0.2*tropoScatterLoss_dB);
    const double val2 = std::pow(10.0, -0.2*modifiedDiffractionAndAnomolousPropagationLoss_dB);
    return -5.0 * std::log10(val1+val2)+ tx_clutterLoss_dB + rx_clutterLoss_dB;
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
