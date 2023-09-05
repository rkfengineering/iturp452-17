#include "MainModel/TropoScatter.h"
#include "MainModel/CalculationHelpers.h"
#include "Common/MathHelpers.h"
#include <filesystem>
#include <iostream>
#include <ostream>
#include <sstream>

double ITUR_P452::TropoScatter::calcTroposcatterLoss_dB(const CommonInputs& commonInputs, const ITUR_P452::TxRxPair&elevationAngles_mrad, 
        const double& eff_radius_med_km, const double& seaLevelSurfaceRefractivity, const double& txHorizonGain_dBi, 
        const double& rxHorizonGain_dBi, const double& temp_K, const double& dryPressure_hPa){

	if (commonInputs.p_percent < 0.001 || commonInputs.p_percent > 50.0) {
        std::ostringstream oStrStream;
        oStrStream << "ERROR: TropoScatter::calcTroposcatterLoss_dB(): " << 
                "The given time percentage falls outside of the range of percentages exceeded ([0.001%-50%]) specified in the ITU-R P452-17 Section 4.3: " 
                << commonInputs.p_percent << " %!" << std::endl;
        throw std::domain_error(oStrStream.str());
	}

    const double path_angular_distance_mrad = Helpers::calcPathAngularDistance_mrad(elevationAngles_mrad,
                commonInputs.d_tot_km,eff_radius_med_km);

    //Equation 45a     
    const double frequencyDependentLoss_dB = 25.0*std::log10(commonInputs.freq_GHz)
                -2.5*ItuModels::MathHelpers::simpleSquare(std::log10(commonInputs.freq_GHz/2.0));
    //Equation 45b
    const double aperatureToMedium_CouplingLoss_dB = 0.051*std::exp(0.055*(txHorizonGain_dBi+rxHorizonGain_dBi));

    //The extra path length from considering the antenna heights is insignificant 
    //but it is also an explicit difference between 452-16 and 452-17
    //Equation 8a line of sight distance accounting for height differential (calculated for transhorizon paths too)
    const double d_los_km = std::sqrt(ItuModels::MathHelpers::simpleSquare(commonInputs.d_tot_km)+
                ItuModels::MathHelpers::simpleSquare((commonInputs.height_tx_asl_m-commonInputs.height_rx_asl_m)/1000.0));
    //Equation 9 using rho = 3 g/m^3
    //The extra path length from considering the antenna heights is insignificant
    const double gasAtten_dB = Helpers::calcGasAtten_dB(d_los_km,commonInputs.freq_GHz,temp_K,dryPressure_hPa,3.0);

    //Equation 45 Empirical troposcatter loss model
    return 190.0 + frequencyDependentLoss_dB + 20.0*std::log10(commonInputs.d_tot_km) + 0.573*path_angular_distance_mrad
            -0.15*seaLevelSurfaceRefractivity + aperatureToMedium_CouplingLoss_dB + gasAtten_dB
            -10.1*std::pow(-std::log10(commonInputs.p_percent/50.0),0.7);
}
