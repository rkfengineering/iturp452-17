#include "MainModel/BasicProp.h"
#include "MainModel/CalculationHelpers.h"
#include "Common/MathHelpers.h"
#include <cmath>

ITUR_P452::BasicProp::BasicProp(const double& d_tot_km, const double& height_tx_asl_m, const double& height_rx_asl_m,
        const double& freq_GHz, const double& temp_K, const double& dryPressure_hPa, const double& frac_over_sea,
        const double& p_percent, const double& b0_percent, const ITUR_P452::TxRxPair& horizonDists_km):
        m_d_tot_km{d_tot_km}, m_height_tx_asl_m{height_tx_asl_m}, m_height_rx_asl_m{height_rx_asl_m},
        m_freq_GHz{freq_GHz}, m_temp_K{temp_K}, m_dryPressure_hPa{dryPressure_hPa},
        m_p_percent{p_percent}, m_b0_percent{b0_percent}, m_horizonDists_km{horizonDists_km},
        m_frac_over_sea{frac_over_sea} {       
}                    

void ITUR_P452::BasicProp::calcTransmissionlosses_dB(double& out_freeSpaceWithGasLoss_dB, double& out_basicTransmissionLoss_p_percent_dB, 
                                double& out_basicTransmissionLoss_b0_percent_dB) const{
    out_freeSpaceWithGasLoss_dB = calcPathLossWithGas_dB();
    //Eq 11
    out_basicTransmissionLoss_p_percent_dB = out_freeSpaceWithGasLoss_dB + calcMultipathFocusingCorrection_dB(m_p_percent);
    //Eq 12
    out_basicTransmissionLoss_b0_percent_dB = out_freeSpaceWithGasLoss_dB + calcMultipathFocusingCorrection_dB(m_b0_percent);    
}

//TODO this code is reused in anomalous prop with the same inputs. consider removing redundant calculations
double ITUR_P452::BasicProp::calcPathLossWithGas_dB() const{
    
    //Equation 8a distance accounting for height differential (calculated for transhorizon paths too)
    const double d_los_km = std::sqrt(MathHelpers::simpleSquare(m_d_tot_km)+
                                        MathHelpers::simpleSquare((m_height_tx_asl_m-m_height_rx_asl_m)/1000.0));
    //Equation 9a Water Vapor Density
    const double rho = 7.5 + 2.5 * m_frac_over_sea;
   
    // Equation 9
    const double gasLoss_dB = Helpers::calcGasAtten_dB(d_los_km,m_freq_GHz,m_temp_K,m_dryPressure_hPa,rho);

    //Equation 8 add free space path loss and gaseous attenuation
    return BasicProp::calcFreeSpacePathLoss_dB(d_los_km,m_freq_GHz) + gasLoss_dB;
}

double ITUR_P452::BasicProp::calcFreeSpacePathLoss_dB(const double& d_los_km, const double& freq_GHz){
    //Equation 8 without gas atten. Modified to avoid calling log10 twice
    return 92.4 + 20.0*std::log10(freq_GHz*d_los_km);
}

double ITUR_P452::BasicProp::calcMultipathFocusingCorrection_dB(const double& p_percent) const{
    //Eq 10a, 10b 
    return 2.6 * (1.0 - std::exp(-0.1*(m_horizonDists_km.first+m_horizonDists_km.second)))*std::log10(p_percent/50.0);
}
