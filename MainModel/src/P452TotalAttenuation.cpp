#include "MainModel/P452TotalAttenuation.h"
#include "MainModel/CalculationHelpers.h"
#include <tuple>

ITUR_P452::TotalClearAirAttenuation::TotalClearAirAttenuation(const double& freq_GHz, const double& p_percent, 
            const PathProfile::Path path_TxToRx, const double& height_tx_m, const double& height_rx_m, 
            const double& centerLatitude_deg, const double& deltaN,
            const ClutterModel::ClutterType& tx_clutterType, const ClutterModel::ClutterType& rx_clutterType):
            m_freq_GHz{freq_GHz}, m_p_percent{p_percent}, m_deltaN{deltaN}{
                pre_calcPathParameters(path_TxToRx,centerLatitude_deg,height_tx_m,height_rx_m,tx_clutterType,rx_clutterType);
}

void ITUR_P452::TotalClearAirAttenuation::pre_calcPathParameters(const PathProfile::Path& path_TxToRx, 
        const double& centerLatitude_deg,const double& height_tx_m, const double& height_rx_m, 
        const ClutterModel::ClutterType& tx_clutterType, const ClutterModel::ClutterType& rx_clutterType){

    //Path Parameters calculated using actual path
    //The conversion from deltaN to median effective radius is performed again in the diffraction model (redundant)
    //The choice is to keep diffraction loss as a standalone module (since it has its own validation data)
    m_effEarthRadius_med_km = Helpers::calcMedianEffectiveRadius_km(m_deltaN);
    m_fracOverSea = path_TxToRx.calcFracOverSea();
    m_b0_percent = path_TxToRx.calcTimePercentBeta0(centerLatitude_deg);

    //Apply height gain model correction from clutter model
    const auto ClutterResults = ClutterModel::calculateClutterModel(m_freq_GHz,path_TxToRx,height_tx_m,height_rx_m,
                                                                    tx_clutterType,rx_clutterType);

    m_mod_path = ClutterResults.modifiedPath;
    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterResults.modifiedHeights_m;
    m_clutterLoss_dB = ClutterResults.clutterLoss_dB;

    m_height_tx_asl_m = hg_height_tx_m + m_mod_path.front().h_asl_m;
    m_height_rx_asl_m = hg_height_rx_m + m_mod_path.back().h_asl_m;
    m_d_tot_km = m_mod_path.back().d_km;

    //Path geometry parameters of modified path
    m_HorizonVals = Helpers::calcHorizonAnglesAndDistances(
        m_mod_path, m_height_tx_asl_m, m_height_rx_asl_m, m_effEarthRadius_med_km, m_freq_GHz
    );
}

double ITUR_P452::TotalClearAirAttenuation::calcTotalClearAirAttenuation(const double& temp_K, const double& dryPressure_hPa, 
    const double& dist_coast_tx_km, const double& dist_coast_rx_km, 
    const double& seaLevelSurfaceRefractivity, const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi,
    const ItuModels::Enumerations::PolarizationType& pol) const{
  
    //submodel calculations--------------------------------------------------------------------------------------
    const auto [HorizonAngles_mrad, HorizonDistances_km] = m_HorizonVals;
    
    double freeSpaceWithGasLoss_dB, basicTransmissionLoss_p_percent_dB, basicTransmissionLoss_b0_percent_dB; 
    const auto BasicPropModel = BasicProp(m_d_tot_km, m_height_tx_asl_m, m_height_rx_asl_m, m_freq_GHz, temp_K, dryPressure_hPa, 
        m_fracOverSea, m_p_percent, m_b0_percent, HorizonDistances_km);
    //Equation 8 (Lbfsg) basic transmission loss with gas atten
    //Equation 11 (Lb0p) basic transmission loss with gas and multipath not exceeded for p percent of time
    //Equation 12 (Lb0b) basic transmission loss with gas and multipath not exceeded for b0 percent of time
    BasicPropModel.calcTransmissionlosses_dB(freeSpaceWithGasLoss_dB, basicTransmissionLoss_p_percent_dB, 
                                            basicTransmissionLoss_b0_percent_dB);

    double diffractionLoss_p_percent_dB, diffractionLoss_median_dB;           
    const auto DiffractionModel = DiffractionLoss(m_mod_path, m_height_tx_asl_m, m_height_rx_asl_m, m_freq_GHz, 
        m_deltaN, pol, m_p_percent, m_b0_percent, m_fracOverSea);
    //Delta Bullington Diffraction Loss calculations
    DiffractionModel.calcDiffractionLoss_dB(diffractionLoss_median_dB,diffractionLoss_p_percent_dB);

    //Anomalous Propagation Calculations (Ducting and Layer Reflection)
    const auto AnomalousPropModel = ITUR_P452::AnomalousProp(m_mod_path, m_freq_GHz, m_height_tx_asl_m, 
        m_height_rx_asl_m, temp_K, dryPressure_hPa, dist_coast_tx_km, dist_coast_rx_km, m_p_percent,
        m_b0_percent, m_effEarthRadius_med_km, m_HorizonVals, m_fracOverSea);
    const double anomalousPropagationLoss_dB = AnomalousPropModel.calcAnomalousPropLoss_dB();
    
    //Calculate Tropospheric Scatter
    const double tropoScatterLoss_dB = TropoScatter::calcTroposcatterLoss_dB(m_d_tot_km,m_freq_GHz,m_height_tx_asl_m,
        m_height_rx_asl_m, HorizonAngles_mrad, m_effEarthRadius_med_km, seaLevelSurfaceRefractivity, 
        txHorizonGain_dBi, rxHorizonGain_dBi, temp_K, dryPressure_hPa, m_p_percent);


    //Total Attenuanuation Calculations------------------------------------------------------------------------------

    //Equation 43 (Lbd50) basic loss with diffraction loss not exceeded for 50 percent of time
    const double basicWithMedianDiffractionLoss_dB = freeSpaceWithGasLoss_dB + diffractionLoss_median_dB;
    //Equation 44 (Lbd) basic loss with diffraction loss not exceeded for p percent of time
    const double basicWithDiffractionLoss_p_percent_dB = basicTransmissionLoss_p_percent_dB + diffractionLoss_p_percent_dB;

    //Equation 60 Minimum basic transmission loss associated with LOS propagation and over-sea sub-path diffraction
    double minLossWithOverSeaSubPathDiffraction_dB = basicTransmissionLoss_p_percent_dB + (1-m_fracOverSea)*diffractionLoss_p_percent_dB;
    if(m_p_percent>=m_b0_percent){
        //Diffraction Interpolation parameter
        const double diffractionInterpolationParameter = 
                    CalculationHelpers::inv_cum_norm(m_p_percent/100.0)/CalculationHelpers::inv_cum_norm(m_b0_percent/100.0);
                    
        minLossWithOverSeaSubPathDiffraction_dB = ItuModels::MathHelpers::interpolate1D(
                basicWithMedianDiffractionLoss_dB, 
                basicTransmissionLoss_b0_percent_dB + (1-m_fracOverSea)*diffractionLoss_p_percent_dB,
                diffractionInterpolationParameter
        );
    }

    //Equation 61 (Lminbap)
    constexpr double eta = 2.5; //constant parameter
    const double minLossWithAnomalousPropagation_dB = 
                eta*std::log(std::exp(anomalousPropagationLoss_dB/eta)+std::exp(basicTransmissionLoss_p_percent_dB/eta));

    //Fk
    const double pathBlendingInterpolationParameter = TotalClearAirAttenuation::calcPathBlendingInterpolationParameter(m_d_tot_km);
    //Equation 62 (Lbda)
    double diffractionAndAnomalousPropagationLoss_dB = basicWithDiffractionLoss_p_percent_dB;
    if(minLossWithAnomalousPropagation_dB <= basicWithDiffractionLoss_p_percent_dB){
        diffractionAndAnomalousPropagationLoss_dB = ItuModels::MathHelpers::interpolate1D(minLossWithAnomalousPropagation_dB,
                                                    basicWithDiffractionLoss_p_percent_dB,pathBlendingInterpolationParameter);
    }

    //Fj
    const double slopeInterpolationParameter = 
        TotalClearAirAttenuation::calcSlopeInterpolationParameter(m_mod_path,m_effEarthRadius_med_km,m_height_tx_asl_m,m_height_rx_asl_m);
    //Equation 63 (Lbam)
    const double modifiedDiffractionAndAnomalousPropagationLoss_dB = 
                                                    ItuModels::MathHelpers::interpolate1D(diffractionAndAnomalousPropagationLoss_dB, 
                                                    minLossWithOverSeaSubPathDiffraction_dB,slopeInterpolationParameter);

    //Equation 64 Total Loss predicted by model, combines losses using a geometric mean of the linear values
    const double val1 = std::pow(10.0, -0.2*tropoScatterLoss_dB);
    const double val2 = std::pow(10.0, -0.2*modifiedDiffractionAndAnomalousPropagationLoss_dB);
    return -5.0 * std::log10(val1+val2)+ m_clutterLoss_dB.first + m_clutterLoss_dB.second;
}

double ITUR_P452::TotalClearAirAttenuation::calcSlopeInterpolationParameter(const PathProfile::Path& path, const double& effEarthRadius_med_km,
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

double ITUR_P452::TotalClearAirAttenuation::calcPathBlendingInterpolationParameter(const double& d_tot_km){
    //fixed parameter for distance range of associated blending
    constexpr double dsw = 20;
    //fixed parameter for blending slope at ends of range
    constexpr double kappa = 0.5;

    //Equation 59 Calculate interpolation factor Fk
    return 1.0 - 0.5 * (1.0 + std::tanh(3.0 * kappa * (d_tot_km-dsw)/dsw));
}
