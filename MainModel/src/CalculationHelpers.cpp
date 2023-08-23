#include "MainModel/CalculationHelpers.h"
#include "Common/PhysicalConstants.h"
#include <cmath>

double ITUR_P452::CalculationHelpers::inv_cum_norm(double prob){
    //if(prob>0.5){
    //    //Invalid input
    //    return std::nan;
    //}

    //min prob value 1e-6
    double tx = std::sqrt(-2*std::log(std::max(prob,1e-6)));
    constexpr double C0 = 2.515516698;
    constexpr double C1 = 0.802853;
    constexpr double C2 = 0.010328;
    constexpr double D1 = 1.432788;
    constexpr double D2 = 0.189269;
    constexpr double D3 = 0.001308;
    
    double ksi = ((C2*tx+C1)*tx+C0)/(((D3*tx+D2)*tx+D1)*tx+1);
    return ksi-tx;
}

double ITUR_P452::CalculationHelpers::convert_freqGHz_to_wavelength_m(const double& freq_GHz){
    return 1e-9*PhysicalConstants::SPEED_OF_LIGHT_M_PER_S/freq_GHz;
}
