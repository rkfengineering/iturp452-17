#include "DiffractionLoss.h"
#include "InvCumNorm.h"
#include "MathHelpers.h"
#include <algorithm>
#include <cmath>

double constexpr SPEED_OF_LIGHT = 0.2998; //m*GHz

double DiffractionLoss::bullLoss(const std::vector<double> d, const std::vector<double> h, const double hts,
        const double hrs, const double ap, const double freqGHz){
    
    const double Ce = 1.0/ap; //effective Earth Curvature
    const double lam = SPEED_OF_LIGHT/freqGHz; //wavelength in m

    double dtot = d.back()-d.front();//total path length
    double Luc = 0;//knife edge loss

    //Find intermediate profile point with highest slope from Tx
    int points = d.size();//assume arrays have intermediate points
    std::vector<double> slopes(points-2);
    for(int i=1; i<points-1; i++){ //TODO optimize later
        slopes[i-1] = (h[i]+500*Ce*d[i]*(dtot-d[i])-hts)/d[i]; //Eq 14
    }
    double Stim = *std::max_element(slopes.begin(),slopes.end());
    //Slope of line from Tx to Rx assuming LOS
    double Str = (hrs-hts)/dtot; //Eq 15

    //Case 1 LOS path
    if(Stim<Str){
        //calculate diffraction parameter at every profile point
        std::vector<double> nus(points-2);
        for(int i=1; i<points-1; i++){ //TODO optimize later
            nus[i-1] = ((h[i]+500*Ce*d[i]*(dtot-d[i])-(hts*(dtot-d[i])+hrs*d[i])/dtot) *
                std::sqrt(0.002*dtot/(lam*d[i]*(dtot-d[i])))); //Eq 16
        }
        double numax = *std::max_element(nus.begin(),nus.end());
  
        if (numax > -0.78){
            //Eq 13, 17 Knife Edge Loss Approximation
            Luc = 6.9 + 20*std::log10(std::sqrt(MathHelpers::simpleSquare(numax-0.1)+1)+numax-0.1);
        }
    }
    //Case 2 Transhorizon path
    else{
        //Find intermediate profile point with highest slope from Rx
        std::vector<double> slopes(points-2);
        for(int i=1; i<points-1; i++){ //TODO optimize later
            slopes[i-1] = (h[i]+500*Ce*d[i]*(dtot-d[i])-hrs)/(dtot-d[i]); //Eq 18
        }
        double Srim = *std::max_element(slopes.begin(),slopes.end());

        double dbp = (hrs-hts+Srim*dtot)/(Stim+Srim); //Eq 19
        double nub = (hts+Stim*dbp-(hts*(dtot-dbp)+hrs*dbp)/dtot) *
            std::sqrt(0.002*dtot/(lam*dbp*(dtot-dbp))); //Eq 20
        
        if (nub > -0.78){
            //Eq 13, 21 Knife Edge Loss Approximation
            Luc = 6.9 + 20*std::log10(std::sqrt(MathHelpers::simpleSquare(nub-0.1)+1)+nub-0.1);
        }

    }

    //Eq 22 Bullington Loss 
    return Luc + (1-std::exp(-Luc/6.0))*(10+0.02*dtot); 
}

double DiffractionLoss::delta_bullington(const ProfilePath& path, const double hts, const double hrs, const double hstd, 
        const double hsrd, const double ap, const double freqGHz,
        const double omega, const Enumerations::PolarizationType pol){

    //TODO figure out if its worth using ProfilePath or just pass two vectors around

    //Actual heights and profile
    double Lbulla = DiffractionLoss::bullLoss(path.d, path.h, hts, hrs, ap, freqGHz);
    
    //modified heights and zero profile
    double hts1 = hts - hstd;
    double hrs1 = hrs - hsrd;
    std::vector<double> h1(path.length,0.0);
    double Lbulls = DiffractionLoss::bullLoss(path.d, h1, hts1, hrs1, ap, freqGHz);

    //Spherical Earth Diffraction Loss
    double dtot = path.d.back() - path.d.front();
    double Ldsph = DiffractionLoss::se_diffLoss(dtot, hts1, hrs1, ap, freqGHz,omega,pol);

    //Eq 40
    return Lbulla + std::max(Ldsph - Lbulls, 0.0);
}

DiffractionLoss::DiffResults DiffractionLoss::diffLoss(const ProfilePath& path, const double hts, const double hrs, const double hstd, 
        const double hsrd, double freqGHz, const double omega, const double p, const double b0, 
        const double DN, const Enumerations::PolarizationType pol){

    //Median value of effective earth radius
    double k50 = 157.0/(157-DN); //Eq 5
    double ae = 6371 * k50; //Eq 6a
    //kb = 3 at Point Incidence of anomalous propagation b0
    double ab = 6371 * 3; //Eq 6b
    
    //Delta Bullington Loss
    double Ld50 = DiffractionLoss::delta_bullington(path,hts,hrs,hstd,hsrd,ae,freqGHz,omega,pol);

    DiffractionLoss::DiffResults results = DiffractionLoss::DiffResults();
    results.Ld50 = Ld50;
    if(p==50){
        results.Ldp = Ld50;
    }
    else if(p<50){
        //Delta Bullington Loss not exceeded for b0% time
        double Ldb = DiffractionLoss::delta_bullington(path,hts,hrs,hstd,hsrd,ab,freqGHz,omega,pol);

        //interpolation factor 
        double Fi;
        if(p>b0){
            Fi = inv_cum_norm(p/100)/inv_cum_norm(b0/100); //Eq 41a
        }
        else{
            Fi = 1;
        }
        results.Ldp = Ld50 + Fi*(Ldb-Ld50);
    }
    //WARNING p>50 is not defined or checked for
    return results;
}

double DiffractionLoss::se_diffLoss(const double d_gc, const double hte, const double hre, const double ap,
        const double freqGHz, const double omega, const Enumerations::PolarizationType pol){

    const double lam = SPEED_OF_LIGHT/freqGHz; //wavelength in m
    //marginal LOS distance for a smooth path
    double dlos = std::sqrt(2*ap)*std::sqrt(0.001*hte) + std::sqrt(0.001*hre);//Eq 23

    //use 4.2.2.1 if applicable
    if(d_gc>=dlos){
        return DiffractionLoss::se_first_term(d_gc,hte,hre,ap,freqGHz,omega,pol);
    }
    //otherwise:

    double c = (hte-hre)/(hte+hre); //Eq 25d
    double m = 250*d_gc*d_gc/(ap*(hte+hre)); //Eq 25e
    double b = 2*std::sqrt((m+1)/(3*m))*std::cos(M_PI/3+
        std::acos(3*c/2*std::sqrt(3* m / MathHelpers::simpleCube(m+1)))/3); //Eq 25c
    double dse1 = d_gc/2*(1+b); //Eq 25a
    double dse2 = d_gc - dse1; //Eq 25b

    double hse = ((hte - 500*dse1*dse1/ap)*dse2 + (hre - 500*dse1*dse2/ap)*dse1)/d_gc; //Eq 24

    //Required Clearance for zero diffraction loss
    double hreq = 17.456 * std::sqrt(dse1*dse2*lam/d_gc); //Eq 26
    if(hse>hreq){
        return 0.0;
    }

    //modified effective earth radius
    double aem = 500*MathHelpers::simpleSquare(d_gc/(std::sqrt(hte)+std::sqrt(hre))); //Eq 27
    //Use 4.2.2.1 method with modified effective earth radius
    double Ldft = DiffractionLoss::se_first_term(d_gc,hte,hre,aem,freqGHz,omega,pol);
    if(Ldft<0){
        return 0.0;
    }

    //Calculate spherical loss by interpolation
    return (1-hse/hreq)*Ldft; //Eq 28
}

double DiffractionLoss::se_first_term_inner(const double epsr, const double sigma, const double d_gc, const double hte, 
        const double hre, const double adft, const double freqGHz, const Enumerations::PolarizationType pol){
    
    //Normalized factor for surface admittance for Horizontal Polarization
    double K = 0.036*std::pow((adft*freqGHz),-1.0/3)*std::pow((MathHelpers::simpleSquare(epsr-1)+
        MathHelpers::simpleSquare(18*sigma/freqGHz)),-1.0/4); //Eq 30a

    //Normalized factor for surface admittance for Vertical Polarization
    if(pol!=Enumerations::PolarizationType::HorizontalPolarized){
        double K_ver = K*std::sqrt(MathHelpers::simpleSquare(epsr)+MathHelpers::simpleSquare(18*sigma/freqGHz)); //Eq 30b
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

    if (freqGHz<0.3){
        double K2 = K*K;
        double K4 = K2*K2;
        beta_dft = (1+1.6*K2 + 0.67*K4)/(1+4.5*K2 + 1.53*K4); //Eq 31
    }

    //Normalized Distance
    double X = 21.88 * beta_dft * std::pow((freqGHz/(adft*adft)),1.0/3) * d_gc; //Eq 32
    
    //Eq 33,36
    double B = 0.9575 * beta_dft*beta_dft * std::pow(freqGHz*freqGHz/adft,1.0/3);
    double Bt = B*hte;
    double Br = B*hre; 

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
    double min_GY = 2+20*std::log10(K);
    GYt = std::max(GYt, min_GY);
    GYr = std::max(GYr, min_GY);

    return -Fx-GYt-GYr; //Eq 37
}         

double DiffractionLoss::se_first_term(const double d_gc, const double hte, const double hre, const double adft,
        const double freqGHz, const double omega, const Enumerations::PolarizationType pol){
   
    //Loss over land, espr = 22, sigma = 0.003
    double Ldft_land = DiffractionLoss::se_first_term_inner(22,0.003,d_gc,hte,hre,adft,freqGHz,pol);

    //Loss over sea, espr = 80, sigma = 5
    double Ldft_sea = DiffractionLoss::se_first_term_inner(80,5,d_gc,hte,hre,adft,freqGHz,pol);

    return omega*Ldft_sea + (1-omega)*Ldft_land; //Eq 29
}