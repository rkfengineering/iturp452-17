#include "ClearAirModel/P452TotalAttenuation.h"
#include "ClearAirModel/CalculationHelpers.h"
#include <tuple>

ClearAirModel::p452_TotalAttenuation::p452_TotalAttenuation(const double& freq_GHz, const double& p_percent, const PathProfile::Path path, 
            const double& height_tx_m, const double& height_rx_m, const double& centerLatitude_deg, const double& txHorizonGain_dBi, 
            const double& rxHorizonGain_dBi, const Enumerations::PolarizationType& pol, const double& dist_coast_tx_km, 
            const double& dist_coast_rx_km, const double& deltaN, const double& surfaceRefractivity,
            const double& temp_K, const double& dryPressure_hPa, const ClutterType& tx_clutterType, 
            const ClutterType& rx_clutterType):
            m_freq_GHz{freq_GHz}, m_p_percent{p_percent}{
                populatePathParameters(path,deltaN,centerLatitude_deg,height_tx_m,height_rx_m,tx_clutterType,rx_clutterType);
                calculateSubModels(temp_K,dryPressure_hPa,deltaN,dist_coast_tx_km,dist_coast_rx_km,
                    surfaceRefractivity,txHorizonGain_dBi,rxHorizonGain_dBi,pol);
                m_totalTransmissionLoss_dB=calcP452TotalAttenuation();
}

void ClearAirModel::p452_TotalAttenuation::populatePathParameters(const PathProfile::Path& path, const double& deltaN, 
        const double& centerLatitude_deg,const double& height_tx_m, const double& height_rx_m, 
        const ClutterType& tx_clutterType, const ClutterType& rx_clutterType){

    //Path Parameters
    m_effEarthRadius_med_km = ClearAirModelHelpers::calcMedianEffectiveRadius_km(deltaN);
    m_fracOverSea = path.calcFracOverSea();
    m_b0_percent = path.calcTimePercentBeta0(centerLatitude_deg);

    //Apply height gain model correction from clutter model
    const auto ClutterLossModel = ClutterLoss(m_freq_GHz,path,height_tx_m,height_rx_m,tx_clutterType,rx_clutterType);
    m_mod_path = ClutterLossModel.getModifiedPath();

    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterLossModel.getHeightGainModelHeights_m();
    std::tie(m_tx_clutterLoss_dB, m_rx_clutterLoss_dB) = ClutterLossModel.getClutterLoss_dB();

    m_height_tx_asl_m = hg_height_tx_m + m_mod_path.front().h_asl_m;
    m_height_rx_asl_m = hg_height_rx_m + m_mod_path.back().h_asl_m;
    m_d_tot_km = m_mod_path.back().d_km;

    //Path geometry parameters of modified path
    m_HorizonVals = ClearAirModelHelpers::calcHorizonAnglesAndDistances(
        m_mod_path, m_height_tx_asl_m, m_height_rx_asl_m, m_effEarthRadius_med_km, m_freq_GHz
    );
}

//TODO replace DN with median effective earth radius as input for diffraction model
void ClearAirModel::p452_TotalAttenuation::calculateSubModels(const double& temp_K, const double& dryPressure_hPa, 
    const double deltaN, const double& dist_coast_tx_km, const double& dist_coast_rx_km, 
    const double& seaLevelSurfaceRefractivity, const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi,
    const Enumerations::PolarizationType& pol){

    const auto [HorizonAngles_mrad, HorizonDistances_km] = m_HorizonVals;
    
    const auto BasicPropModel = BasicProp(m_d_tot_km, m_height_tx_asl_m, m_height_rx_asl_m, m_freq_GHz, temp_K, dryPressure_hPa, 
        m_fracOverSea, m_p_percent, m_b0_percent, HorizonDistances_km);
    //Equation 8 (Lbfsg)
    m_freeSpaceWithGasLoss_dB = BasicPropModel.getFreeSpaceWithGasLoss_dB();
    //Equation 11 (Lb0p) basic transmission loss with gas and multipath not exceeded for p percent of time
    m_basicTransmissionLoss_p_percent_dB = BasicPropModel.getBasicTransmissionLoss_p_percent_dB();
    //Equation 12 (Lb0b)
    m_basicTransmissionLoss_b0_percent_dB = BasicPropModel.getBasicTransmissionLoss_b0_percent_dB();

    const auto DiffractionModel = DiffractionLoss(m_mod_path, m_height_tx_asl_m, m_height_rx_asl_m, m_freq_GHz, 
        deltaN, pol, m_p_percent, m_b0_percent, m_fracOverSea);
    //Delta Bullington Diffraction Loss calculations
    m_diffractionLoss_p_percent_dB = DiffractionModel.getDiffractionLoss_p_percent_dB();
    m_diffractionLoss_median_dB = DiffractionModel.getDiffractionLoss_median_dB();

    //Anomolous Propagation Calculations (Ducting and Layer Reflection)
    const auto AnomolousPropModel = ClearAirModel::AnomolousProp(m_mod_path, m_freq_GHz, m_height_tx_asl_m, 
        m_height_rx_asl_m, temp_K, dryPressure_hPa, dist_coast_tx_km, dist_coast_rx_km, m_p_percent,
        m_b0_percent, m_effEarthRadius_med_km, m_HorizonVals, m_fracOverSea);
    m_anomolousPropagationLoss_dB = AnomolousPropModel.getAnomolousPropLoss_dB();
    
    //Calculate Tropospheric Scatter
    m_tropoScatterLoss_dB = TropoScatter::calcTroposcatterLoss_dB(m_d_tot_km,m_freq_GHz,m_height_tx_asl_m,
        m_height_rx_asl_m, HorizonAngles_mrad, m_effEarthRadius_med_km, seaLevelSurfaceRefractivity, 
        txHorizonGain_dBi, rxHorizonGain_dBi, temp_K, dryPressure_hPa, m_p_percent);
}

double ClearAirModel::p452_TotalAttenuation::calcP452TotalAttenuation(){
  
    //Equation 43 (Lbd50) basic loss with diffraction loss not exceeded for 50 percent of time
    const double basicWithMedianDiffractionLoss_dB = m_freeSpaceWithGasLoss_dB + m_diffractionLoss_median_dB;
    //Equation 44 (Lbd) basic loss with diffraction loss not exceeded for p percent of time
    const double basicWithDiffractionLoss_p_percent_dB = m_basicTransmissionLoss_p_percent_dB + m_diffractionLoss_p_percent_dB;

    //Equation 60 Minimum basic transmission loss associated with LOS propagation and over-sea sub-path diffraction
    double minLossWithOverSeaSubPathDiffraction_dB = m_basicTransmissionLoss_p_percent_dB + (1-m_fracOverSea)*m_diffractionLoss_p_percent_dB;
    if(m_p_percent>=m_b0_percent){
        //Diffraction Interpolation parameter
        const double diffractionInterpolationParameter = 
                    CalculationHelpers::inv_cum_norm(m_p_percent/100.0)/CalculationHelpers::inv_cum_norm(m_b0_percent/100.0);
                    
        minLossWithOverSeaSubPathDiffraction_dB = MathHelpers::interpolate1D(
                basicWithMedianDiffractionLoss_dB, 
                m_basicTransmissionLoss_b0_percent_dB + (1-m_fracOverSea)*m_diffractionLoss_p_percent_dB,
                diffractionInterpolationParameter
        );
    }

    //Equation 61 (Lminbap)
    constexpr double eta = 2.5; //constant parameter
    const double minLossWithAnomolousPropagation_dB = 
                eta*std::log(std::exp(m_anomolousPropagationLoss_dB/eta)+std::exp(m_basicTransmissionLoss_p_percent_dB/eta));

    //Fk
    const double pathBlendingInterpolationParameter = p452_TotalAttenuation::calcPathBlendingInterpolationParameter(m_d_tot_km);
    //Equation 62 (Lbda)
    double diffractionAndAnomolousPropagationLoss_dB = basicWithDiffractionLoss_p_percent_dB;
    if(minLossWithAnomolousPropagation_dB <= basicWithDiffractionLoss_p_percent_dB){
        diffractionAndAnomolousPropagationLoss_dB = MathHelpers::interpolate1D(minLossWithAnomolousPropagation_dB,
                                                    basicWithDiffractionLoss_p_percent_dB,pathBlendingInterpolationParameter);
    }

    //Fj
    const double slopeInterpolationParameter = 
        p452_TotalAttenuation::calcSlopeInterpolationParameter(m_mod_path,m_effEarthRadius_med_km,m_height_tx_asl_m,m_height_rx_asl_m);
    //Equation 63 (Lbam)
    const double modifiedDiffractionAndAnomolousPropagationLoss_dB = 
                                                    MathHelpers::interpolate1D(diffractionAndAnomolousPropagationLoss_dB, 
                                                    minLossWithOverSeaSubPathDiffraction_dB,slopeInterpolationParameter);

    //Equation 64 Total Loss predicted by model, combines losses using a geometric mean of the linear values
    const double val1 = std::pow(10.0, -0.2*m_tropoScatterLoss_dB);
    const double val2 = std::pow(10.0, -0.2*modifiedDiffractionAndAnomolousPropagationLoss_dB);
    return -5.0 * std::log10(val1+val2)+ m_tx_clutterLoss_dB + m_rx_clutterLoss_dB;
}

double ClearAirModel::p452_TotalAttenuation::calcSlopeInterpolationParameter(const PathProfile::Path path, const double effEarthRadius_med_km,
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

double ClearAirModel::p452_TotalAttenuation::calcPathBlendingInterpolationParameter(const double& d_tot_km){
    //fixed parameter for distance range of associated blending
    constexpr double dsw = 20;
    //fixed parameter for blending slope at ends of range
    constexpr double kappa = 0.5;

    //Equation 59 Calculate interpolation factor Fk
    return 1.0 - 0.5 * (1.0 + std::tanh(3.0 * kappa * (d_tot_km-dsw)/dsw));
}
