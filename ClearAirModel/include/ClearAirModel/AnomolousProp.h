#ifndef ANOMOLOUS_PROP_H
#define ANOMOLOUS_PROP_H

#include "PathProfile.h"
#include "EffectiveEarth.h"

//Section 4.4 Prediction of the basic transmission loss, Lba (dB) 
//occurring during periods of anomalous propagation (ducting and layer reflection)
/*
namespace AnomolousProp {
    //TODO try and break this up more
    double anomolousPropLoss(const double& d_tot_km, const EffectiveEarth::HorizonAnglesAndDistances& horizonVals,
            const double& dist_coast_tx_km, const double& dist_coast_rx_km, const double& height_tx_m_amsl, 
            const double& height_rx_m_amsl, const EffectiveEarth::TxRxPair& effHeights_ducting_m, const double& terrainRoughness_m,
            const double& freqGHz, const double& p_percent, const double& temp_K, const double& dryPressure_hPa, const double& frac_over_sea,
            const double& eff_radius_med_km, const double& b0);

}//end namespace AnomolousProp*/

#endif /* ANOMOLOUS_PROP_H */


/*
    
    //Equation 8a distance accounting for height differential
    const double d_los_km = std::sqrt(MathHelpers::simpleSquare(d_tot_km)+
                                        MathHelpers::simpleSquare((height_tx_asl_m-height_rx_asl_m)/1000.0));
    //Equation 9a Water Vapor Density
    const double rho = 7.5 + 2.5 * frac_over_sea;
   
    // Equation 9
    const double gasLoss_dB = BasicProp::gasAttenWrapper_dB(d_los_km,freq_GHz,temp_K,dryPressure_hPa,rho);
    
    */
