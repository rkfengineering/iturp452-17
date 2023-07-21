#include "ClearAirModel/TropoScatter.h"
#include "Common/MathHelpers.h"
#include <filesystem>
#include <iostream>
#include <ostream>
#include <sstream>

//TODO find a good place to put this
//data grid object to store surface refractivity (N0) data
static const DataGrid2 seaLevelSurfaceRefractivityMap{CMAKE_CLEARAIR_SRC_DIR / std::filesystem::path("data/N050.TXT"),1.5};

double TropoScatter::calcTroposcatterLoss_dB(const double& d_tot_km, const double& freq_GHz, const double& path_angular_distance_mrad, 
        const double& seaLevelSurfaceRefractivity, const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi, 
        const double& temp_K, const double& dryPressure_hPa, const double& p_percent){

	if (p_percent < 0.001 || p_percent > 50.0) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: TropoScatter::calcTroposcatterLoss_dB(): " << 
			"The given time percentage falls outside of the range of percentages exceeded ([0.001%-50%]) specified in the ITU-R P452-17 Section 4.3: " 
			<< p_percent << " %!" << std::endl;
		throw std::domain_error(oStrStream.str());
	}

    //Equation 45a     
    const double frequencyDependentLoss_dB = 25.0*std::log10(freq_GHz)-2.5*MathHelpers::simpleSquare(std::log10(freq_GHz/2.0));
    //Equation 45b
    const double aperatureToMedium_CouplingLoss_dB = 0.051*std::exp(0.055*(txHorizonGain_dBi+rxHorizonGain_dBi));

    //Equation 9 using rho = 3 g/m^3
    //The extra path length from considering the antenna heights is insignificant
    const double gasAtten_dB = BasicProp::calcGasAtten_dB(d_tot_km,freq_GHz,temp_K,dryPressure_hPa,3.0);

    //Equation 45 Empirical troposcatter loss model
    return 190.0 + frequencyDependentLoss_dB + 20.0*std::log10(d_tot_km) + 0.573*path_angular_distance_mrad
            -0.15*seaLevelSurfaceRefractivity + aperatureToMedium_CouplingLoss_dB + gasAtten_dB
            -10.1*std::pow(-std::log10(p_percent/50.0),0.7);
}

double TropoScatter::fetchSeaLevelSurfaceRefractivity(const GeodeticCoord& location){
    return seaLevelSurfaceRefractivityMap.interpolate2D(location);
}