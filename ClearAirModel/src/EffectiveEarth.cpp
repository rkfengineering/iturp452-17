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

EffectiveEarth::TxRxPair EffectiveEarth::calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m(const PathProfile::Path& path, 
        const double& height_tx_m, const double& height_rx_m){

    const double d_tot = path.back().d_km; //assume distances start at 0
    const double height_tx_m_amsl = path.front().h_asl_m + height_tx_m;
    const double height_rx_m_amsl = path.back().h_asl_m + height_rx_m;

    //Section 5.1.6.3

    //Highest Obstruction height (m)
    double height_obs_max = std::numeric_limits<double>::lowest();
    //Horizon elevation angles (mrad)
    double alpha_obs_t_max = std::numeric_limits<double>::lowest();
    double alpha_obs_r_max = std::numeric_limits<double>::lowest();

    //get max values for intermediate obstruction
    double height_val,delta_d;
    PathProfile::ProfilePoint point;
    for (auto cit = path.begin()+1; cit<path.end()-1; ++cit){
        point = *cit;
        delta_d = d_tot-point.d_km;
        height_val = point.h_asl_m-(height_tx_m_amsl*delta_d+height_rx_m_amsl*point.d_km)/d_tot; //Eq 165d
        height_obs_max = std::max(height_obs_max, height_val);//Eq 165a
        alpha_obs_t_max = std::max(alpha_obs_t_max,height_val/point.d_km);             //Eq 165b
        alpha_obs_r_max = std::max(alpha_obs_r_max,height_val/delta_d);      //Eq 165c
    }

    //Equations 166a, 166b
    //Tx,Rx heights from a least squares smooth path
    auto [height_smooth_tx_amsl_m,height_smooth_rx_amsl_m] = 
            EffectiveEarth::calcLeastSquaresSmoothEarthTxRxHeights_helper_amsl_m(path);

    //Modify heights to compensate for obstructions
    if(height_obs_max>0){
        const double v1 = alpha_obs_t_max+alpha_obs_r_max;
        const double gt = alpha_obs_t_max/v1; //Eq 166e
        const double gr = alpha_obs_r_max/v1; //Eq 166f

        height_smooth_tx_amsl_m-=height_obs_max*gt; //Eq 166c
        height_smooth_rx_amsl_m-=height_obs_max*gr; //Eq 166d
    }

    //Limit effective antenna heights to be above actual terrain ground height
    double eff_height_tx_amsl_m = std::min(path.front().h_asl_m, height_smooth_tx_amsl_m); //Eq 167 a,b
    double eff_height_rx_amsl_m = std::min(path.back().h_asl_m, height_smooth_rx_amsl_m); //Eq 167 c,d
    
    return EffectiveEarth::TxRxPair{eff_height_tx_amsl_m,eff_height_rx_amsl_m};
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

EffectiveEarth::TxRxPair EffectiveEarth::calcSmoothEarthTxRxHeights_DuctingModel_amsl_m(const PathProfile::Path& path, 
        const double& height_tx_m, const double& height_rx_m){

    //Equations 166a, 166b
    //Tx,Rx heights from a least squares smooth path
    auto [height_smooth_tx_amsl_m,height_smooth_rx_amsl_m] = 
            EffectiveEarth::calcLeastSquaresSmoothEarthTxRxHeights_helper_amsl_m(path);
    //Equation 168 terminal heights must be above ground level
    height_smooth_tx_amsl_m = std::min(height_smooth_tx_amsl_m, path.front().h_asl_m);
    height_smooth_rx_amsl_m = std::min(height_smooth_rx_amsl_m, path.back().h_asl_m);

    //Equation 170
    const double eff_height_tx_m = height_tx_m + path.front().h_asl_m - height_smooth_tx_amsl_m;
    const double eff_height_rx_m = height_rx_m + path.back().h_asl_m - height_smooth_rx_amsl_m;

    return EffectiveEarth::TxRxPair(eff_height_tx_m,eff_height_rx_m);
}

//TODO refactor code to make more efficient. We don't have to call the helper function this often
double EffectiveEarth::calcTerrainRoughness_m(const PathProfile::Path& path, TxRxPair horizonDists_km){
    //Equations 166a, 166b
    //Tx,Rx heights from a least squares smooth path
    auto [height_smooth_tx_amsl_m,height_smooth_rx_amsl_m] = 
            EffectiveEarth::calcLeastSquaresSmoothEarthTxRxHeights_helper_amsl_m(path);
    //Equation 168 terminal heights must be above ground level
    height_smooth_tx_amsl_m = std::min(height_smooth_tx_amsl_m, path.front().h_asl_m);
    height_smooth_rx_amsl_m = std::min(height_smooth_rx_amsl_m, path.back().h_asl_m);

    //smooth earth surface slope
    //assume path starts at 0 km 
    const double slope = (height_smooth_rx_amsl_m-height_smooth_tx_amsl_m)/path.back().d_km;

    //only evaluate section between horizon points
    const auto [tx_horizon_km, rx_horizon_from_rx] = horizonDists_km;
    const double rx_horizon_km = path.back().d_km -rx_horizon_from_rx;

    //Equation 171 calculate terrain roughness above smooth earth path
    double terrainRoughness_m = 0; //the parameter can never be negative 
    double heightAboveSmoothPath;
    for(auto point : path){
        if(point.d_km>=tx_horizon_km && point.d_km<=rx_horizon_km){
            heightAboveSmoothPath = point.h_asl_m-(height_smooth_tx_amsl_m + slope*point.d_km);
            terrainRoughness_m = std::max(terrainRoughness_m, heightAboveSmoothPath);
        }
    }
    return terrainRoughness_m;
}


