#include "EffectiveEarth.h"
#include "PhysicalConstants.h"
#include <limits>
#include <cmath>

double EffectiveEarth::eff_radius_p50_km(const double& delta_N){
    const double k50 = 157.0/(157.0-delta_N); //Eq 5
    return 6371.0 * k50; //Eq 6a
}

EffectiveEarth::TxRxPair EffectiveEarth::smoothEarthHeights_AMSL(const PathProfile::Path& path){

    const double d_tot = path.back().d_km; //assume distances start at 0
    //Section 5.1.6.2
    double v1 = 0;
    double v2 = 0;
    PathProfile::ProfilePoint prevPoint = path.front();
    for(uint32_t i = 1; i<path.size(); ++i){ //start at second point
        const PathProfile::ProfilePoint point = path.at(i);
        v1+=(point.d_km-prevPoint.d_km)*(point.h_masl+prevPoint.h_masl);//Eq 161
        v2+=(point.d_km-prevPoint.d_km)*(point.h_masl*(2*point.d_km+prevPoint.d_km)+prevPoint.h_masl*(point.d_km+2*prevPoint.d_km)); //Eq 162
        prevPoint = point;
    }
    
    EffectiveEarth::TxRxPair height_results = EffectiveEarth::TxRxPair();

    height_results.tx_val = (2.0*v1*d_tot-v2)/(d_tot*d_tot); //Eq 163
    height_results.rx_val = (v2-v1*d_tot)/(d_tot*d_tot);   //Eq 164

    return height_results;
}

EffectiveEarth::TxRxPair EffectiveEarth::smoothEarthHeights_diffractionModel(const PathProfile::Path& path, 
        const double& height_tx_m, const double& height_rx_m){

    const double d_tot = path.back().d_km; //assume distances start at 0
    const double height_tx_m_amsl = path.front().h_masl + height_tx_m;
    const double height_rx_m_amsl = path.back().h_masl + height_rx_m;

    //Section 5.1.6.3

    //Highest Obstruction height (m)
    double height_obs_max = std::numeric_limits<double>::lowest();
    //Horizon elevation angles (mrad)
    double alpha_obs_t_max = std::numeric_limits<double>::lowest();
    double alpha_obs_r_max = std::numeric_limits<double>::lowest();

    //get max values
    double height_val,delta_d;
    PathProfile::ProfilePoint point;
    for (auto cit = path.begin()+1; cit<path.end()-1; ++cit){
        point = *cit;
        delta_d = d_tot-point.d_km;
        height_val = point.h_masl-(height_tx_m_amsl*delta_d+height_rx_m_amsl*point.d_km)/d_tot; //Eq 165d
        height_obs_max = std::max(height_obs_max, height_val);//Eq 165a
        alpha_obs_t_max = std::max(alpha_obs_t_max,height_val/point.d_km);             //Eq 165b
        alpha_obs_r_max = std::max(alpha_obs_r_max,height_val/delta_d);      //Eq 165c
    }

    const auto height_se_m_amsl_pair = EffectiveEarth::smoothEarthHeights_AMSL(path);

    double hstp = height_se_m_amsl_pair.tx_val; //Eq 166a
    double hsrp = height_se_m_amsl_pair.rx_val; //Eq 166b
    
    if(height_obs_max>0){
        const double v1 = alpha_obs_t_max+alpha_obs_r_max;
        const double gt = alpha_obs_t_max/v1; //Eq 166e
        const double gr = alpha_obs_r_max/v1; //Eq 166f

        hstp-=height_obs_max*gt; //Eq 166c
        hsrp-=height_obs_max*gr; //Eq 166d
    }
    
    EffectiveEarth::TxRxPair height_results = EffectiveEarth::TxRxPair();
    height_results.tx_val = std::min(path.front().h_masl, hstp); //Eq 167 a,b
    height_results.rx_val = std::min(path.back().h_masl, hsrp); //Eq 167 c,d
    return height_results;
}

EffectiveEarth::HorizonAnglesAndDistances EffectiveEarth::getHorizonAnglesAndDistances(const PathProfile::Path& path,
            const double& height_tx_masl, const double& height_rx_masl, const double& eff_radius_med_km, const double& freq_GHz){

    const double d_tot = path.back().d_km;
    //Equation 153 Angle from tx to rx, relative to local horizon
    const double theta_td = 1e3*std::atan(
        (height_rx_masl-height_tx_masl)/(1e3*d_tot)
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
            (pt.h_masl-height_tx_masl)/(1e3*pt.d_km)
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

    HorizonAnglesAndDistances results = HorizonAnglesAndDistances();

    //Equation 154 Tx (interfering) antenna horizon elevation angle
    results.horizonElevation_mrad.tx_val = std::max(theta_tmax,theta_td);

    if(isTranshorizon){
        //Equation 155
        results.horizonDist_km.tx_val = path.at(tx_index).d_km;

        //Calculate max rx elevation angle
        double theta_rmax = std::numeric_limits<double>::lowest();
        double theta_j,delta_d;
        uint16_t rx_index=0;
        for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
            delta_d = d_tot-cit->d_km;
            //Equation 157 calculate elevation angle from rx to terrain point
            theta_j = 1e3*std::atan(
                (cit->h_masl-height_rx_masl)/(1e3*delta_d)
                -delta_d/(2.0*eff_radius_med_km)
            );
            //assume prefer points closer to rx
            if(theta_j>=theta_rmax){
                theta_rmax = theta_j;
                rx_index = cit - path.cbegin();
            }
        }

        //Equation 156b horizon elevation angle from rx antenna
        results.horizonElevation_mrad.rx_val = theta_rmax;
        //Equation 158
        results.horizonDist_km.rx_val = d_tot - path.at(rx_index).d_km;

    }
    else{
        //find bullington point for LOS path

        //Effective Earth Curvature
        const double Ce = 1.0/eff_radius_med_km;
        //wavelength in m
        const double lam = PhysicalConstants::SPEED_OF_LIGHT_M_GHZ/freq_GHz;

        double v1,v2,nu,delta_d;
        //calculate diffraction parameter at every intermediate profile point 
        double numax = std::numeric_limits<double>::lowest();
        for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
            delta_d = d_tot-cit->d_km;

            //Eq 16,155a function to calculate diffraction parameter nu
            v1 = (cit->h_masl+500.0*Ce*cit->d_km*(delta_d)-(height_tx_masl*(delta_d)+height_rx_masl*cit->d_km)/d_tot);
            v2 = std::sqrt(0.002*d_tot/(lam*cit->d_km*delta_d));
            nu = v1*v2;
            //assume prefer points closer to tx
            if(nu>numax){
                numax = nu;
                tx_index = cit - path.cbegin();
            }
        }
        //Equation 155a
        results.horizonDist_km.tx_val = path.at(tx_index).d_km;

        //Equation 156a
        results.horizonElevation_mrad.rx_val = 1e3*std::atan(
            (height_tx_masl-height_rx_masl)/(1e3*d_tot)
            -d_tot/(2.0*eff_radius_med_km)
        );
        //Equation 158a
        results.horizonDist_km.rx_val = d_tot - results.horizonDist_km.tx_val;
    }
    return results;
}


/*
EffectiveEarth::SmoothEarthResults EffectiveEarth::smoothEarthHeights(const PathProfile::Path& path, double& height_tx_m,
        double& height_rx_m, double& eff_radius_med_km, double& freq_GHz){

    const double d_tot = path.back().d_km; //assume distances start at 0
    const double hts = path.front().h_masl + height_tx_m;
    const double hrs = path.back().h_masl + height_rx_m;

    //Section 5.1.6.2
    double v1 = 0;
    double v2 = 0;
    PathProfile::ProfilePoint prevPoint = path.front();
    for(uint32_t i = 1; i<path.size(); i++){ //start at second point
        const PathProfile::ProfilePoint point = path.at(i);
        v1+=(point.d_km-prevPoint.d_km)*(point.h_masl+prevPoint.h_masl);//Eq 161
        v2+=(point.d_km-prevPoint.d_km)*(point.h_masl*(2*point.d_km-prevPoint.d_km)+prevPoint.h_masl*(point.d_km-2*prevPoint.d_km)); //Eq 162
        prevPoint = point;
    }
    
    const double height_tx_se_m_amsl = (2*v1*d_tot-v2)/(d_tot*d_tot); //Eq 163
    const double height_rx_se_m_amsl = (v2-v1*d_tot)/(d_tot*d_tot);   //Eq 164

    //Section 5.1.6.3
    double HHi = path.h[1]-(hts*(d_tot-path.d[1])+hrs*path.d[1])/d_tot;
    double hobs = HHi;
    double alpha_obt = HHi/path.d[1];
    double alpha_obr = HHi/(d_tot-path.d[1]);

    //get max values
    for(int i=2; i<path.length-1;i++){
        HHi = path.h[i]-(hts*(d_tot-path.d[i])+hrs*path.d[i])/d_tot; //Eq 165d
        hobs = std::max(hobs, HHi);                                //Eq 165a
        alpha_obt = std::max(alpha_obt,HHi/path.d[i]);             //Eq 165b
        alpha_obr = std::max(alpha_obr,HHi/(d_tot-path.d[i]));      //Eq 165c
    }

    double gt = alpha_obt/(alpha_obt+alpha_obr); //Eq 166e
    double gr = alpha_obr/(alpha_obt+alpha_obr); //Eq 166f

    double hstp = height_tx_se_m_amsl; //Eq 166a
    double hsrp = height_rx_se_m_amsl; //Eq 166b
    if(hobs>0){
        hstp-=hobs*gt; //Eq 166c
        hsrp-=hobs*gr; //Eq 166d
    }
    
    SmoothEarthResults results = SmoothEarthResults();
    results.height_tx_se_m_amsl = height_tx_se_m_amsl;
    results.height_rx_se_m_amsl = height_rx_se_m_amsl;
    results.eff_height_tx_diffraction_m = std::min(path.h.front(), hstp); //Eq 167 a,b
    results.eff_height_rx_diffraction_m = std::min(path.h.back(), hsrp); //Eq 167 c,d



    //horizon elevation angle and distance
    double theta_max;
    int kindex=0;
    double theta = 1000*std::atan((path.h[1]-hts)/1000*path.d[1])-path.d[1]/(2*eff_radius_med_km);
    theta_max = theta;
    for(int i=2; i<path.length-1;i++){
        theta = 1000*std::atan((path.h[i]-hts)/1000*path.d[i])-path.d[i]/(2*eff_radius_med_km);
        if(theta>theta_max){//in case of ties, favor the point closer to tx
            theta_max = theta;
            kindex = i;
        }
    }

    for(int i=2; i<path.length-1;i++){
        theta = 1000*std::atan((path.h[i]-hts)/1000*path.d[i])-path.d[i]/(2*eff_radius_med_km);
        if(theta>=theta_max){//in case of ties, favor the point closer to rx
            theta_max = theta;
            kindex = i;
        }
    }

    return results;
}*/