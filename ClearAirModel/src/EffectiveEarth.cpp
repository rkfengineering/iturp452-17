#include "ClearAirModel/EffectiveEarth.h"
#include "ClearAirModel/CalculationHelpers.h"
#include "Common/PhysicalConstants.h"
#include <limits>
#include <cmath>
#include <iostream>
#include <ostream>
#include <sstream>

double EffectiveEarth::calcMedianEffectiveRadius_km(const double& delta_N){
    const double k50 = 157.0/(157.0-delta_N); //Eq 5
    return 6371.0 * k50; //Eq 6a
}

//Least Squares linear approximation of the actual path
EffectiveEarth::TxRxPair EffectiveEarth::calcLeastSquaresSmoothEarthTxRxHeights_helper_amsl_m(const PathProfile::Path& path){

    const double d_tot = path.back().d_km; //assume distances start at 0
    //Section 5.1.6.2
    double v1 = 0;
    double v2 = 0;
    PathProfile::ProfilePoint prevPoint = path.front();
    for(uint32_t i = 1; i<path.size(); ++i){ //start at second point
        const PathProfile::ProfilePoint point = path.at(i);
        //Equation 161
        v1+=(point.d_km-prevPoint.d_km)*(point.h_asl_m+prevPoint.h_asl_m);
        //Equation 162
        v2+=(point.d_km-prevPoint.d_km)*(point.h_asl_m*(2*point.d_km+prevPoint.d_km)+prevPoint.h_asl_m*(point.d_km+2*prevPoint.d_km)); //Eq 162
        prevPoint = point;
    }
    
    //Equation 163 Tx least squares model height 
    const double height_tx_amsl_m = (2.0*v1*d_tot-v2)/(d_tot*d_tot); 
    //Equation 164 Tx least squares model height 
    const double height_rx_amsl_m = (v2-v1*d_tot)/(d_tot*d_tot);   

    return EffectiveEarth::TxRxPair{height_tx_amsl_m,height_rx_amsl_m};
}

EffectiveEarth::HorizonAnglesAndDistances EffectiveEarth::calcHorizonAnglesAndDistances(const PathProfile::Path& path,
            const double& height_tx_asl_m, const double& height_rx_asl_m, const double& eff_radius_med_km, const double& freq_GHz){

    const double d_tot = path.back().d_km;
    //Equation 153 Angle from tx to rx, relative to local horizon
    const double theta_td = 1e3*std::atan(
        (height_rx_asl_m-height_tx_asl_m)/(1e3*d_tot)
        -d_tot/(2.0*eff_radius_med_km)
    );

    // calculate tx elevation angle
    double theta_tmax = std::numeric_limits<double>::lowest();
    
    uint16_t tx_index=0;
    double theta_i;
    PathProfile::ProfilePoint pt;
    //Eq 151 max elevation angle from tx to terrain point
    for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
        pt = *cit;
        //Equation 152 function to calculate elevation angle from tx to terrain point
        theta_i = 1e3*std::atan(
            (pt.h_asl_m-height_tx_asl_m)/(1e3*pt.d_km)
            -pt.d_km/(2.0*eff_radius_med_km)
        );
        //assume prefer points closer to tx
        if(theta_i>theta_tmax){
            theta_tmax = theta_i;
            tx_index = cit - path.cbegin();
        }
    }
    
    //Equation 150 check if path is Line of Sight or Trans-Horizon
    const bool isTranshorizon = theta_tmax>theta_td;

    //declare variables for final results
    double horizonElevation_tx_mrad, horizonElevation_rx_mrad, horizonDist_tx_km, horizonDist_rx_km;

    //Equation 154 Tx (interfering) antenna horizon elevation angle
    horizonElevation_tx_mrad = std::max(theta_tmax,theta_td);

    if(isTranshorizon){
        //Equation 155
        horizonDist_tx_km = path.at(tx_index).d_km;

        //Calculate max rx elevation angle
        double theta_rmax = std::numeric_limits<double>::lowest();
        double theta_j,delta_d;
        uint16_t rx_index=0;
        for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
            delta_d = d_tot-cit->d_km;
            //Equation 157 calculate elevation angle from rx to terrain point
            theta_j = 1e3*std::atan(
                (cit->h_asl_m-height_rx_asl_m)/(1e3*delta_d)
                -delta_d/(2.0*eff_radius_med_km)
            );
            //assume prefer points closer to rx
            if(theta_j>=theta_rmax){
                theta_rmax = theta_j;
                rx_index = cit - path.cbegin();
            }
        }

        //Equation 156b horizon elevation angle from rx antenna
        horizonElevation_rx_mrad = theta_rmax;
        //Equation 158
        horizonDist_rx_km = d_tot - path.at(rx_index).d_km;

    }
    else{
        //find bullington point for LOS path

        //Effective Earth Curvature
        const double Ce = 1.0/eff_radius_med_km;
        //wavelength in m
        const double wavelength_m =  CalculationHelpers::convert_freqGHz_to_wavelength_m(freq_GHz);

        double v1,v2,nu,delta_d;
        //calculate diffraction parameter at every intermediate profile point 
        double numax = std::numeric_limits<double>::lowest();
        for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
            delta_d = d_tot-cit->d_km;

            //Eq 16,155a function to calculate diffraction parameter nu
            v1 = (cit->h_asl_m+500.0*Ce*cit->d_km*(delta_d)-(height_tx_asl_m*(delta_d)+height_rx_asl_m*cit->d_km)/d_tot);
            v2 = std::sqrt(0.002*d_tot/(wavelength_m*cit->d_km*delta_d));
            nu = v1*v2;
            //assume prefer points closer to tx
            if(nu>numax){
                numax = nu;
                tx_index = cit - path.cbegin();
            }
        }
        //Equation 155a
        horizonDist_tx_km = path.at(tx_index).d_km;

        //Equation 156a
        horizonElevation_rx_mrad = 1e3*std::atan(
            (height_tx_asl_m-height_rx_asl_m)/(1e3*d_tot)
            -d_tot/(2.0*eff_radius_med_km)
        );
        //Equation 158a
        horizonDist_rx_km = d_tot - horizonDist_tx_km;
    }

    return EffectiveEarth::HorizonAnglesAndDistances{
        EffectiveEarth::TxRxPair{horizonElevation_tx_mrad,horizonElevation_rx_mrad},
        EffectiveEarth::TxRxPair{horizonDist_tx_km,horizonDist_rx_km}
    };
}

double EffectiveEarth::calcPathAngularDistance_mrad(const EffectiveEarth::TxRxPair& elevationAngles_mrad, 
        const double& dtot_km, const double& eff_radius_med_km){

    const auto [txHorizonAngle_mrad, rxHorizonAngle_mrad] = elevationAngles_mrad;
    //Equation 159
    return 1e3*dtot_km/eff_radius_med_km + txHorizonAngle_mrad + rxHorizonAngle_mrad;
}


