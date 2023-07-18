#include "DiffractionLoss.h"
#include "InvCumNorm.h"
#include "MathHelpers.h"
#include "EffectiveEarth.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include "PhysicalConstants.h"

double DiffractionLoss::bullLoss(const PathProfile::Path& path, const double& height_tx_masl,
        const double& height_rx_masl, const double& eff_radius_p_km, const double& freq_GHz){
    
    const double Ce = 1.0/eff_radius_p_km; //effective Earth Curvature
    const double lam = PhysicalConstants::SPEED_OF_LIGHT_M_GHZ/freq_GHz; //wavelength in m

    const double d_tot = path.back().d_km-path.front().d_km;//total path length
    double loss_knifeEdge_dB = 0;//knife edge loss

    //Find intermediate profile point with highest slope from Tx

    //need to exclude first and last point

    //Eq 14 get max slope to profile point from tx 
    double max_slope_tx = std::numeric_limits<double>::lowest();
    double slope_tx;
    PathProfile::ProfilePoint pt_tx;
    for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
        pt_tx = *cit;
        slope_tx = (pt_tx.h_masl+500*Ce*pt_tx.d_km*(d_tot-pt_tx.d_km)-height_tx_masl)/pt_tx.d_km;
        max_slope_tx = std::max(max_slope_tx,slope_tx); 
    }

    //Eq 15 Slope of line from Tx to Rx assuming LOS
    const double slope_tr_los = (height_rx_masl-height_tx_masl)/d_tot;

    //Case 1 LOS path
    if(max_slope_tx<slope_tr_los){

        //calculate diffraction parameter at every intermediate profile point 
        double numax = std::numeric_limits<double>::lowest();
        double v1,v2,delta_d;
        PathProfile::ProfilePoint pt;

        //Eq 16 function to calculate diffraction parameter nu
        //Also see Eq 155a
        for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
            pt = *cit;
            delta_d = d_tot-pt.d_km;
            //Note. There is no floor operation in this equation. The square brackets may not be rendered correctly. See Eq 155a
            v1 = (pt.h_masl+500.0*Ce*pt.d_km*(delta_d)-(height_tx_masl*(delta_d)+height_rx_masl*pt.d_km)/d_tot);
            v2 = std::sqrt(0.002*d_tot/(lam*pt.d_km*delta_d));
            numax = std::max(numax,v1*v2); 
        }

        if (numax > -0.78){
            //Eq 13, 17 Knife Edge Loss Approximation
            loss_knifeEdge_dB = 6.9 + 20.0*std::log10(std::sqrt(MathHelpers::simpleSquare(numax-0.1)+1.0)+numax-0.1);
        }
    }
    //Case 2 Transhorizon path
    else{
        //Find intermediate profile point with highest slope from Rx

        //Eq 18 get max slope to profile point from rx 
        double max_slope_rx = std::numeric_limits<double>::lowest();
        double slope_rx;
        PathProfile::ProfilePoint pt_rx;
        for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
            pt_rx = *cit;
            slope_rx = (pt_rx.h_masl+500*Ce*pt_rx.d_km*(d_tot-pt_rx.d_km)-height_rx_masl)/(d_tot-pt_rx.d_km); 
            max_slope_rx = std::max(max_slope_rx,slope_rx); 
        }

        //Eq 19 distance from bullington point to tx
        const double dbp = (height_rx_masl-height_tx_masl+max_slope_rx*d_tot)/(max_slope_tx+max_slope_rx); 

        //Eq 20 diffraction parameter nu at bullington point
        const double nub = (height_tx_masl+max_slope_tx*dbp-(height_tx_masl*(d_tot-dbp)+height_rx_masl*dbp)/d_tot) *
            std::sqrt(0.002*d_tot/(lam*dbp*(d_tot-dbp)));
        
        if (nub > -0.78){
            //Eq 13, 21 Knife Edge Loss Approximation
            loss_knifeEdge_dB = 6.9 + 20.0*std::log10(std::sqrt(MathHelpers::simpleSquare(nub-0.1)+1.0)+nub-0.1);
        }
    }
    //Eq 22 Bullington Loss 
    return loss_knifeEdge_dB + (1-std::exp(-loss_knifeEdge_dB/6.0))*(10+0.02*d_tot); 
}

double DiffractionLoss::delta_bullington(const PathProfile::Path& path, const double& height_tx_masl, const double& height_rx_masl,
        const double& eff_height_itx_masl, const double& eff_height_irx_masl, const double& eff_radius_p_km, const double& freq_GHz,
        const double& frac_over_sea, const Enumerations::PolarizationType& pol){

    //Actual heights and profile
    const double Lbulla = DiffractionLoss::bullLoss(path, height_tx_masl, height_rx_masl, eff_radius_p_km, freq_GHz);
    
    //modified heights and zero profile
    PathProfile::Path zeroHeightPath;
    for (auto point : path){
        zeroHeightPath.push_back(PathProfile::ProfilePoint(point.d_km, 0.0));
    }
    const double mod_height_tx_masl= height_tx_masl- eff_height_itx_masl;
    const double mod_height_rx_masl= height_rx_masl- eff_height_irx_masl;
    const double Lbulls = DiffractionLoss::bullLoss(zeroHeightPath, mod_height_tx_masl, mod_height_rx_masl, eff_radius_p_km, freq_GHz);

    //Spherical Earth Diffraction Loss
    const double d_tot = path.back().d_km - path.front().d_km;
    const double Ldsph = DiffractionLoss::se_diffLoss(d_tot, mod_height_tx_masl, mod_height_rx_masl, 
                                                eff_radius_p_km, freq_GHz,frac_over_sea,pol);

    //Eq 40 Delta Bullington Diffraction Loss (dB)
    return Lbulla + std::max(Ldsph - Lbulls, 0.0);
}

DiffractionLoss::DiffResults DiffractionLoss::diffLoss(const PathProfile::Path& path, const double& height_tx_masl, 
        const double& height_rx_masl, const double& eff_height_itx_masl, const double& eff_height_irx_masl, 
        const double& freq_GHz, const double& frac_over_sea, const double& p_percent, const double& b0, 
        const double& DN, const Enumerations::PolarizationType& pol){

    const double val_eff_radius_p50_km = EffectiveEarth::eff_radius_p50_km(DN);
    
    //Delta Bullington Loss
    const double diff_loss_p50_dB = DiffractionLoss::delta_bullington(path,height_tx_masl,height_rx_masl,
            eff_height_itx_masl,eff_height_irx_masl,val_eff_radius_p50_km,freq_GHz,frac_over_sea,pol);

    DiffractionLoss::DiffResults results = DiffractionLoss::DiffResults();
    results.diff_loss_p50_dB = diff_loss_p50_dB;
    if(p_percent==50){
        results.diff_loss_p_dB = diff_loss_p50_dB;
    }
    else if(p_percent<50){
        //Delta Bullington Loss not exceeded for b0% time
        const double Ldb = DiffractionLoss::delta_bullington(path,height_tx_masl,height_rx_masl,eff_height_itx_masl,
                eff_height_irx_masl,EffectiveEarth::eff_radius_pb_km,freq_GHz,frac_over_sea,pol);

        //interpolation factor 
        double Fi;
        if(p_percent>b0){
            Fi = inv_cum_norm(p_percent/100)/inv_cum_norm(b0/100); //Eq 41a
        }
        else{
            Fi = 1;
        }
        results.diff_loss_p_dB = diff_loss_p50_dB + Fi*(Ldb-diff_loss_p50_dB);
    }
    //WARNING p>50 is not defined or checked for
    return results;
}

double DiffractionLoss::se_diffLoss(const double& distance_gc_km, const double& eff_height_itx_m, 
        const double& eff_height_irx_m, const double& eff_radius_p_km,
        const double& freq_GHz, const double& frac_over_sea, const Enumerations::PolarizationType& pol){

    const double lam = PhysicalConstants::SPEED_OF_LIGHT_M_GHZ/freq_GHz; //wavelength in m
    //Equation 23 marginal LOS distance for a smooth path
    const double d_los_km = std::sqrt(2.0*eff_radius_p_km)*(std::sqrt(0.001*eff_height_itx_m) + std::sqrt(0.001*eff_height_irx_m));

    //use 4.2.2.1 if applicable
    if(distance_gc_km>=d_los_km){
        return DiffractionLoss::se_first_term(distance_gc_km,eff_height_itx_m,eff_height_irx_m,eff_radius_p_km,freq_GHz,frac_over_sea,pol);
    }
    //otherwise:

    const double c = (eff_height_itx_m-eff_height_irx_m)/(eff_height_itx_m+eff_height_irx_m); //Eq 25d
    const double m = 250*distance_gc_km*distance_gc_km/(eff_radius_p_km*(eff_height_itx_m+eff_height_irx_m)); //Eq 25e
    const double b = 2*std::sqrt((m+1)/(3*m))*std::cos(M_PI/3+
        std::acos(3*c/2*std::sqrt(3* m / MathHelpers::simpleCube(m+1)))/3); //Eq 25c
    const double dse1 = distance_gc_km/2*(1+b); //Eq 25a
    const double dse2 = distance_gc_km - dse1; //Eq 25b

    const double h_se = ((eff_height_itx_m - 500.0*dse1*dse1/eff_radius_p_km)*dse2 + 
                        (eff_height_irx_m - 500.0*dse2*dse2/eff_radius_p_km)*dse1)
                        /distance_gc_km; //Eq 24

    //Equation 26 Required Clearance for zero diffraction loss
    const double h_req_m = 17.456 * std::sqrt(dse1*dse2*lam/distance_gc_km);
    if(h_se>h_req_m){
        return 0.0;
    }

    //modified effective earth radius
    const double aem = 500*MathHelpers::simpleSquare(distance_gc_km/(std::sqrt(eff_height_itx_m)+std::sqrt(eff_height_irx_m))); //Eq 27
    //Use 4.2.2.1 method with modified effective earth radius
    const double loss_firstTerm_dB = DiffractionLoss::se_first_term(distance_gc_km,eff_height_itx_m,
                                                                    eff_height_irx_m,aem,freq_GHz,frac_over_sea,pol);
    
    if(loss_firstTerm_dB<0.0){
        return 0.0;
    }

    //Calculate spherical loss by interpolation (dB)
    return (1.0-h_se/h_req_m)*loss_firstTerm_dB; //Eq 28
}

double DiffractionLoss::se_first_term_inner(const double& eps_r, const double& sigma, const double& distance_gc_km, 
        const double& eff_height_itx_m, const double& eff_height_irx_m, const double& eff_radius_km,
        const double& freq_GHz, const Enumerations::PolarizationType& pol){
    
    //Normalized factor for surface admittance for Horizontal Polarization
    double K = 0.036*std::pow((eff_radius_km*freq_GHz),-1.0/3.0)*std::pow((MathHelpers::simpleSquare(eps_r-1.0)+
        MathHelpers::simpleSquare(18.0*sigma/freq_GHz)),-1.0/4.0); //Eq 30a

    //Normalized factor for surface admittance for Vertical Polarization
    if(pol!=Enumerations::PolarizationType::HorizontalPolarized){
        const double K_ver = K*std::sqrt(MathHelpers::simpleSquare(eps_r)+MathHelpers::simpleSquare(18*sigma/freq_GHz)); //Eq 30b
        if(pol == Enumerations::PolarizationType::VerticalPolarized){
            K = K_ver;
        }
        else{
            //Circular Polarization
            //decompose into components and combine results by vector sum of field amplitude
            //not tested, edge case unlikely to be used
            K = std::sqrt(K*K+K_ver*K_ver);
        }
    } 

    //We can use beta_dft=1 for frequencies above 300 MHz
    double beta_dft=1;

    if (freq_GHz<0.3){
        const double K2 = K*K;
        const double K4 = K2*K2;
        beta_dft = (1.0+1.6*K2 + 0.67*K4)/(1.0+4.5*K2 + 1.53*K4); //Eq 31
    }

    //Normalized Distance
    const double X = 21.88 * beta_dft * std::pow((freq_GHz/(eff_radius_km*eff_radius_km)),1.0/3) * distance_gc_km; //Eq 32
    
    //Eq 33,36
    const double Y = 0.9575 *beta_dft * std::pow(freq_GHz*freq_GHz/eff_radius_km,1.0/3);
    const double Yt = Y*eff_height_itx_m;
    const double Yr = Y*eff_height_irx_m;
    const double Bt = beta_dft*Yt;//This can be optimized for speed if needed
    const double Br = beta_dft*Yr;//Yt and Yr are exposed for easier debugging

    //Distance Term, Eq 34
    double Fx;
    if(X>=1.6){
        Fx = 11+10*std::log10(X)-17.6*X;
    }
    else{
        Fx = -20*std::log10(X) - 5.6488*std::pow(X,1.425);
    }

    //Normalized Height Function, Eq 35
    double GYt, GYr;
    if(Bt>2){
        GYt = 17.6*std::sqrt(Bt-1.1)-5*std::log10(Bt-1.1)-8;
    }
    else{
        GYt = 20*std::log10(Bt+0.1*MathHelpers::simpleCube(Bt));
    }
    
    if(Br>2){
        GYr = 17.6*std::sqrt(Br-1.1)-5*std::log10(Br-1.1)-8;
    }
    else{
        GYr = 20*std::log10(Br+0.1*MathHelpers::simpleCube(Br));
    } 

    //enforce minimum values
    const double min_GY = 2+20*std::log10(K);
    GYt = std::max(GYt, min_GY);
    GYr = std::max(GYr, min_GY);

    return -Fx-GYt-GYr; //Eq 37
}         

double DiffractionLoss::se_first_term(const double& distance_gc_km, const double& eff_height_itx_m, 
        const double& eff_height_irx_m, const double& eff_radius_km,
        const double& freq_GHz, const double& frac_over_sea, const Enumerations::PolarizationType& pol){
   
    //Loss over land, espr = 22, sigma = 0.003
    const double loss_firstTerm_land_dB = DiffractionLoss::se_first_term_inner(22,0.003,distance_gc_km,eff_height_itx_m,
                                                                    eff_height_irx_m, eff_radius_km,freq_GHz,pol);

    //Loss over sea, espr = 80, sigma = 5
    const double loss_firstTerm_sea_dB = DiffractionLoss::se_first_term_inner(80,5,distance_gc_km,eff_height_itx_m,
                                                                    eff_height_irx_m,eff_radius_km,freq_GHz,pol);

    return frac_over_sea*loss_firstTerm_sea_dB + (1-frac_over_sea)*loss_firstTerm_land_dB; //Eq 29
}