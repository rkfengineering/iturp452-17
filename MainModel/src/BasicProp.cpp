#include "MainModel/BasicProp.h"
#include "MainModel/CalculationHelpers.h"
#include "Common/MathHelpers.h"
#include <cmath>

ITUR_P452::BasicProp::BasicProp(const CommonInputs& commonInputs, const double& temp_K, const double& dryPressure_hPa, 
        const ITUR_P452::TxRxPair& horizonDists_km):
        m_commonInputs{commonInputs}, m_temp_K{temp_K}, m_dryPressure_hPa{dryPressure_hPa},
        m_horizonDists_km{horizonDists_km}{       
}                    

void ITUR_P452::BasicProp::calcTransmissionlosses_dB(double& out_freeSpaceWithGasLoss_dB, double& out_basicTransmissionLoss_p_percent_dB, 
                                double& out_basicTransmissionLoss_b0_percent_dB) const{
    out_freeSpaceWithGasLoss_dB = calcPathLossWithGas_dB();
    //Eq 11
    out_basicTransmissionLoss_p_percent_dB = out_freeSpaceWithGasLoss_dB + calcMultipathFocusingCorrection_dB(m_commonInputs.p_percent);
    //Eq 12
    out_basicTransmissionLoss_b0_percent_dB = out_freeSpaceWithGasLoss_dB + 
                calcMultipathFocusingCorrection_dB(m_commonInputs.timePercentBeta0);    
}

//TODO this code is reused in anomalous prop with the same inputs. consider removing redundant calculations
double ITUR_P452::BasicProp::calcPathLossWithGas_dB() const{
    
    //Equation 8a distance accounting for height differential (calculated for transhorizon paths too)
    const double d_los_km = std::sqrt(ItuModels::MathHelpers::simpleSquare(m_commonInputs.d_tot_km)+
                ItuModels::MathHelpers::simpleSquare((m_commonInputs.height_tx_asl_m-m_commonInputs.height_rx_asl_m)/1000.0));
    //Equation 9a Water Vapor Density
    const double rho = 7.5 + 2.5 * m_commonInputs.fracOverSea;
   
    // Equation 9
    const double gasLoss_dB = Helpers::calcGasAtten_dB(d_los_km,m_commonInputs.freq_GHz,m_temp_K,m_dryPressure_hPa,rho);

    //Equation 8 add free space path loss and gaseous attenuation
    return BasicProp::calcFreeSpacePathLoss_dB(d_los_km,m_commonInputs.freq_GHz) + gasLoss_dB;
}

double ITUR_P452::BasicProp::calcFreeSpacePathLoss_dB(const double& d_los_km, const double& freq_GHz){
    //Equation 8 without gas atten. Modified to avoid calling log10 twice
    return 92.4 + 20.0*std::log10(freq_GHz*d_los_km);
}

double ITUR_P452::BasicProp::calcMultipathFocusingCorrection_dB(const double& p_percent) const{
    //Eq 10a, 10b 
    return 2.6 * (1.0 - std::exp(-0.1*(m_horizonDists_km.first+m_horizonDists_km.second)))*std::log10(p_percent/50.0);
}
