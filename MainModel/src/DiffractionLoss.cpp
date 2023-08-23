#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <ostream>
#include <sstream>
#include "Common/PhysicalConstants.h"
#include "Common/MathHelpers.h"
#include "MainModel/DiffractionLoss.h"
#include "MainModel/CalculationHelpers.h"
#include "MainModel/Helpers.h"

ITUR_P452::DiffractionLoss::DiffractionLoss(const PathProfile::Path& path, const double& height_tx_asl_m, const double& height_rx_asl_m,
    const double& freq_GHz, const double& deltaN, const Enumerations::PolarizationType& pol, 
    const double& p_percent, const double&b0_percent, const double& frac_over_sea):
    m_path{path}, m_height_tx_asl_m{height_tx_asl_m}, m_height_rx_asl_m{height_rx_asl_m},
    m_freq_GHz{freq_GHz}, m_deltaN{deltaN}, m_pol{pol},
    m_p_percent{p_percent}, m_b0_percent{b0_percent}, m_frac_over_sea{frac_over_sea} {

    //Path Calculations
    m_d_tot_km = m_path.back().d_km;
    //effective heights for smooth path
    const auto [eff_terrainHeight_itx_asl_m,eff_terrainHeight_irx_asl_m] = calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m();
    m_eff_height_itx_m = m_height_tx_asl_m - eff_terrainHeight_itx_asl_m;
    m_eff_height_irx_m = m_height_rx_asl_m - eff_terrainHeight_irx_asl_m;   
}

void ITUR_P452::DiffractionLoss::calcDiffractionLoss_dB(double& out_diff_loss_median_dB, double& out_diff_loss_p_percent_dB) const{
    const double medianEffectiveRadius_km = Helpers::calcMedianEffectiveRadius_km(m_deltaN);
    out_diff_loss_median_dB = calcDeltaBullingtonLoss_dB(medianEffectiveRadius_km);

    if(m_p_percent>50 || m_p_percent < 0.001){
		std::ostringstream oStrStream;
		oStrStream << "ERROR: DiffractionLoss::calcDiffractionLoss_p_percent_dB(): " << 
			"The given time percentage p falls outside of the range ([0.001-50] %) specified: " 
			<< m_p_percent << " %!" << std::endl;
		throw std::domain_error(oStrStream.str());
	}

    if(m_p_percent==50){
        out_diff_loss_p_percent_dB = out_diff_loss_median_dB;
    }
    else{ //m_p_percent<50
        //Delta Bullington Loss not exceeded for b0_percent% time
        const double diffractionLoss_b0percent_dB = calcDeltaBullingtonLoss_dB(Helpers::k_eff_radius_bpercentExceeded_km);

        //interpolation factor 
        double Fi = 1;
        if(m_p_percent>m_b0_percent){
            Fi = CalculationHelpers::inv_cum_norm(m_p_percent/100.0)/CalculationHelpers::inv_cum_norm(m_b0_percent/100.0); //Eq 41a
        }

        out_diff_loss_p_percent_dB = MathHelpers::interpolate1D(out_diff_loss_median_dB, diffractionLoss_b0percent_dB, Fi);
    }
}

double ITUR_P452::DiffractionLoss::calcDeltaBullingtonLoss_dB(const double& eff_radius_p_km) const{

    //Bullington Loss for the Actual Terrain
    const double Lbulla = calcBullingtonLoss_dB(m_path, m_height_tx_asl_m, m_height_rx_asl_m, eff_radius_p_km);
    
    //modified heights and zero profile
    PathProfile::Path zeroHeightpath;
    for (auto point : m_path){
        zeroHeightpath.push_back(PathProfile::ProfilePoint(point.d_km, 0.0));
    }

    //Bullington Loss for an equivalent Smooth Earth m_path
    const double Lbulls = calcBullingtonLoss_dB(zeroHeightpath, m_eff_height_itx_m, m_eff_height_irx_m, eff_radius_p_km);

    //Spherical Earth Diffraction Loss
    const double Ldsph = calcSphericalEarthDiffractionLoss_dB(eff_radius_p_km);

    //Eq 40 Delta Bullington Diffraction Loss (dB)
    return Lbulla + std::max(Ldsph - Lbulls, 0.0);
}

double ITUR_P452::DiffractionLoss::calcBullingtonLoss_dB(const PathProfile::Path& path, const double& height_tx_asl_m,
        const double& height_rx_asl_m, const double& eff_radius_p_km) const{
    
    const double Ce = 1.0/eff_radius_p_km; //effective Earth Curvature
    const double wavelength_m = CalculationHelpers::convert_freqGHz_to_wavelength_m(m_freq_GHz);
    double loss_knifeEdge_dB = 0;//knife edge loss

    //Find intermediate profile point with highest slope from Tx

    //need to exclude first and last point

    //Eq 14 get max slope to profile point from tx 
    double max_slope_tx = std::numeric_limits<double>::lowest();
    double slope_tx;
    PathProfile::ProfilePoint pt_tx;
    for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
        pt_tx = *cit;
        slope_tx = (pt_tx.h_asl_m+500*Ce*pt_tx.d_km*(m_d_tot_km-pt_tx.d_km)-height_tx_asl_m)/pt_tx.d_km;
        max_slope_tx = std::max(max_slope_tx,slope_tx); 
    }

    //Eq 15 Slope of line from Tx to Rx assuming LOS
    const double slope_tr_los = (height_rx_asl_m-height_tx_asl_m)/m_d_tot_km;
    //Case 1 LOS m_path
    if(max_slope_tx<slope_tr_los){

        //calculate diffraction parameter at every intermediate profile point 
        double numax = std::numeric_limits<double>::lowest();
        double v1,v2,delta_d;
        PathProfile::ProfilePoint pt;

        //Eq 16 function to calculate diffraction parameter nu
        //Also see Eq 155a
        for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
            pt = *cit;
            delta_d = m_d_tot_km-pt.d_km;
            //Note. There is no floor operation in this equation. The square brackets may not be rendered correctly. See Eq 155a
            v1 = (pt.h_asl_m+500.0*Ce*pt.d_km*(delta_d)-(height_tx_asl_m*(delta_d)+height_rx_asl_m*pt.d_km)/m_d_tot_km);
            v2 = std::sqrt(0.002*m_d_tot_km/(wavelength_m*pt.d_km*delta_d));
            numax = std::max(numax,v1*v2); 
        }

        if (numax > -0.78){
            //Eq 13, 17 Knife Edge Loss Approximation
            loss_knifeEdge_dB = 6.9 + 20.0*std::log10(std::sqrt(MathHelpers::simpleSquare(numax-0.1)+1.0)+numax-0.1);
        }
    }
    //Case 2 Transhorizon m_path
    else{
        //Find intermediate profile point with highest slope from Rx

        //Eq 18 get max slope to profile point from rx 
        double max_slope_rx = std::numeric_limits<double>::lowest();
        double slope_rx;
        PathProfile::ProfilePoint pt_rx;
        for(auto cit = path.cbegin()+1; cit<path.cend()-1;++cit){
            pt_rx = *cit;
            slope_rx = (pt_rx.h_asl_m+500*Ce*pt_rx.d_km*(m_d_tot_km-pt_rx.d_km)-height_rx_asl_m)/(m_d_tot_km-pt_rx.d_km); 
            max_slope_rx = std::max(max_slope_rx,slope_rx); 
        }

        //Eq 19 distance from bullington point to tx
        const double dbp = (height_rx_asl_m-height_tx_asl_m+max_slope_rx*m_d_tot_km)/(max_slope_tx+max_slope_rx); 

        //Eq 20 diffraction parameter nu at bullington point
        const double nub = (height_tx_asl_m+max_slope_tx*dbp-(height_tx_asl_m*(m_d_tot_km-dbp)+height_rx_asl_m*dbp)/m_d_tot_km) *
            std::sqrt(0.002*m_d_tot_km/(wavelength_m*dbp*(m_d_tot_km-dbp)));

        if (nub > -0.78){
            //Eq 13, 21 Knife Edge Loss Approximation
            loss_knifeEdge_dB = 6.9 + 20.0*std::log10(std::sqrt(MathHelpers::simpleSquare(nub-0.1)+1.0)+nub-0.1);
        }
    }
    //Eq 22 Bullington Loss 
    return loss_knifeEdge_dB + (1-std::exp(-loss_knifeEdge_dB/6.0))*(10+0.02*m_d_tot_km); 
}

double ITUR_P452::DiffractionLoss::calcSphericalEarthDiffractionLoss_dB(const double& eff_radius_p_km) const{

    const double wavelength_m =  CalculationHelpers::convert_freqGHz_to_wavelength_m(m_freq_GHz); //wavelength in m
    //Equation 23 marginal LOS distance for a smooth m_path
    const double d_los_km = std::sqrt(2.0*eff_radius_p_km)*(std::sqrt(0.001*m_eff_height_itx_m) + std::sqrt(0.001*m_eff_height_irx_m));

    //use 4.2.2.1 if applicable
    if(m_d_tot_km>=d_los_km){
        return calcSphericalEarthDiffraction_firstTerm_dB(eff_radius_p_km);
    }
    //otherwise:

    const double c = (m_eff_height_itx_m-m_eff_height_irx_m)/(m_eff_height_itx_m+m_eff_height_irx_m); //Eq 25d
    const double m = 250*m_d_tot_km*m_d_tot_km/(eff_radius_p_km*(m_eff_height_itx_m+m_eff_height_irx_m)); //Eq 25e
    const double b = 2*std::sqrt((m+1)/(3*m))*std::cos(M_PI/3+
        std::acos(3*c/2*std::sqrt(3* m / MathHelpers::simpleCube(m+1)))/3); //Eq 25c
    const double dse1 = m_d_tot_km/2*(1+b); //Eq 25a
    const double dse2 = m_d_tot_km - dse1; //Eq 25b

    const double h_se = ((m_eff_height_itx_m - 500.0*dse1*dse1/eff_radius_p_km)*dse2 + 
                        (m_eff_height_irx_m - 500.0*dse2*dse2/eff_radius_p_km)*dse1)
                        /m_d_tot_km; //Eq 24

    //Equation 26 Required Clearance for zero diffraction loss
    const double h_req_m = 17.456 * std::sqrt(dse1*dse2*wavelength_m/m_d_tot_km);
    if(h_se>h_req_m){
        return 0.0;
    }

    //modified effective earth radius
    const double mod_effEarthRadius_km = 500*MathHelpers::simpleSquare(m_d_tot_km/(std::sqrt(m_eff_height_itx_m)+std::sqrt(m_eff_height_irx_m))); //Eq 27
    //Use 4.2.2.1 method with modified effective earth radius
    const double loss_firstTerm_dB = calcSphericalEarthDiffraction_firstTerm_dB(mod_effEarthRadius_km);
    
    if(loss_firstTerm_dB<0.0){
        return 0.0;
    }

    //Calculate spherical loss by interpolation (dB)
    return (1.0-h_se/h_req_m)*loss_firstTerm_dB; //Eq 28
}

double ITUR_P452::DiffractionLoss::calcSphericalEarthDiffraction_firstTerm_dB(const double& eff_radius_km) const{

    //Loss over land, relative permittivity = 22, conductivity = 0.003 S/m
    const double loss_firstTerm_land_dB = calcSphericalEarthDiffraction_firstTerm_singleZone_dB(22,0.003,eff_radius_km);

    //Loss over sea, relative permittivity = 80, conductivity = 5 S/m
    const double loss_firstTerm_sea_dB = calcSphericalEarthDiffraction_firstTerm_singleZone_dB(80,5,eff_radius_km);

    //Equation 29
    return MathHelpers::interpolate1D(loss_firstTerm_land_dB, loss_firstTerm_sea_dB, m_frac_over_sea);
}


double ITUR_P452::DiffractionLoss::calcSphericalEarthDiffraction_firstTerm_singleZone_dB(const double& relPermittivity, 
                                                    const double& conductivity,const double& eff_radius_km) const{
    
    //Normalized factor for surface admittance for Horizontal Polarization
    double K = 0.036*std::pow((eff_radius_km*m_freq_GHz),-1.0/3.0)*std::pow((MathHelpers::simpleSquare(relPermittivity-1.0)+
        MathHelpers::simpleSquare(18.0*conductivity/m_freq_GHz)),-1.0/4.0); //Eq 30a

    //Normalized factor for surface admittance for Vertical Polarization
    if(m_pol!=Enumerations::PolarizationType::HorizontalPolarized){
        //Equation 30b
        const double K_ver = K*std::sqrt(
                    MathHelpers::simpleSquare(relPermittivity)
                    + MathHelpers::simpleSquare(18.0*conductivity/m_freq_GHz)
        );
        if(m_pol == Enumerations::PolarizationType::VerticalPolarized){
            K = K_ver;
        }
        else{
            //Circular Polarization
            //decompose into components and combine results by vector sum of field amplitude
            //not tested, edge case unlikely to be used
            K = std::sqrt(K*K+K_ver*K_ver);
            
            std::ostringstream oStrStream;
            std::cerr << "WARNING: DiffractionLoss::calcSphericalEarthDiffraction_firstTerm_singleZone_dB(): " 
            << "This method has not been tested with Circular Polarization";
        }
    } 

    //We can use beta_dft=1 for frequencies above 300 MHz

    const double K2 = K*K;
    const double K4 = K2*K2;
    const double beta_dft = (1.0+1.6*K2 + 0.67*K4)/(1.0+4.5*K2 + 1.53*K4); //Eq 31

    //Normalized Distance
    const double X = 21.88 * beta_dft * std::pow((m_freq_GHz/(eff_radius_km*eff_radius_km)),1.0/3) * m_d_tot_km; //Eq 32
    
    //Eq 33,36
    const double Y = 0.9575 *beta_dft * std::pow(m_freq_GHz*m_freq_GHz/eff_radius_km,1.0/3);
    const double Yt = Y*m_eff_height_itx_m;
    const double Yr = Y*m_eff_height_irx_m;
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

    //Equation 35 Normalized Height Function
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

ITUR_P452::TxRxPair ITUR_P452::DiffractionLoss::calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m() const{

    const double d_tot = m_path.back().d_km; //assume distances start at 0

    //Section 5.1.6.3

    //Highest Obstruction height (m)
    double height_obs_max = std::numeric_limits<double>::lowest();
    //Horizon elevation angles (mrad)
    double alpha_obs_t_max = std::numeric_limits<double>::lowest();
    double alpha_obs_r_max = std::numeric_limits<double>::lowest();

    //get max values for intermediate obstruction
    double height_val,delta_d;
    PathProfile::ProfilePoint point;
    for (auto cit = m_path.begin()+1; cit<m_path.end()-1; ++cit){
        point = *cit;
        delta_d = d_tot-point.d_km;
        height_val = point.h_asl_m-(m_height_tx_asl_m*delta_d+m_height_rx_asl_m*point.d_km)/d_tot; //Eq 165d
        height_obs_max = std::max(height_obs_max, height_val);//Eq 165a
        alpha_obs_t_max = std::max(alpha_obs_t_max,height_val/point.d_km);             //Eq 165b
        alpha_obs_r_max = std::max(alpha_obs_r_max,height_val/delta_d);      //Eq 165c
    }

    //Equations 166a, 166b
    //Tx,Rx heights from a least squares smooth m_path
    auto [height_smooth_tx_amsl_m,height_smooth_rx_amsl_m] = 
            Helpers::calcLeastSquaresSmoothEarthTxRxHeights_helper_amsl_m(m_path);

    //Modify heights to compensate for obstructions
    if(height_obs_max>0){
        const double v1 = alpha_obs_t_max+alpha_obs_r_max;
        const double gt = alpha_obs_t_max/v1; //Eq 166e
        const double gr = alpha_obs_r_max/v1; //Eq 166f

        height_smooth_tx_amsl_m-=height_obs_max*gt; //Eq 166c
        height_smooth_rx_amsl_m-=height_obs_max*gr; //Eq 166d
    }

    //Limit effective antenna heights to be above actual terrain ground height
    double eff_height_tx_amsl_m = std::min(m_path.front().h_asl_m, height_smooth_tx_amsl_m); //Eq 167 a,b
    double eff_height_rx_amsl_m = std::min(m_path.back().h_asl_m, height_smooth_rx_amsl_m); //Eq 167 c,d
    
    return ITUR_P452::TxRxPair{eff_height_tx_amsl_m,eff_height_rx_amsl_m};
}