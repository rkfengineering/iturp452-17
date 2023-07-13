#include "PathProfile.h"
#include <fstream>
#include <sstream>
#include <cmath>

//constructors
PathProfile::ProfilePoint::ProfilePoint(){
}
PathProfile::ProfilePoint::ProfilePoint(double distance_km, double height_masl):
    d_km{distance_km}, h_masl{height_masl}{
}
PathProfile::ProfilePoint::ProfilePoint(double distance_km, double height_masl, ZoneType zonetype):
    d_km{distance_km}, h_masl{height_masl}, zone{zonetype}{
}

//Default Constructor
PathProfile::Path::Path(){
}

//TODO look for a more standard approach to interfacing with csvs
PathProfile::Path::Path(std::string csvPath){
    std::ifstream file;
    file.open(csvPath);
    std::string line;
    std::getline(file,line);//throw away first header line

    std::stringstream ss(line);
    std::string val;
    u_int16_t n = 0;
    while(std::getline(ss,val,',')){
        n+=1;
    }

    while(std::getline(file,line)){
        std::stringstream ss(line);
        std::string val;
        PathProfile::ProfilePoint point;

        for(u_int16_t i=0; i<n; i++){
            std::getline(ss,val,',');
            if(i==0){
                point.d_km = std::stod(val);
            }
            else if(i==1){
                point.h_masl = std::stod(val);
            }
            else if(i==2){
                point.zone = static_cast<PathProfile::ZoneType>(std::stoi(val));
            }
        }
        push_back(point);
    }
}

/*
SmoothEarthResults smoothEarthHeights(const ProfilePath& path, double htg, double hrg, double ae, double freqGHz){
    double dtot = path.d.back(); //assume distances start at 0
    double hts = path.h.front() + htg;
    double hrs = path.h.back() + hrg;

    //Section 5.1.6.2
    double v1 = 0;
    double v2 = 0;
    for(int i = 1; i<path.length; i++){
        v1+=(path.d[i]-path.d[i-1])*(path.h[i]+path.h[i-1]);//Eq 161
        v2+=(path.d[i]-path.d[i-1])*(path.h[i]*(2*path.d[i]+path.d[i-1])+path.h[i-1]*(path.d[i]+2*path.d[i-1])); //Eq 162
    }
    
    double hst = (2*v1*dtot-v2)/(dtot*dtot); //Eq 163
    double hsr = (v2-v1*dtot)/(dtot*dtot);   //Eq 164

    //Section 5.1.6.3
    double HHi = path.h[1]-(hts*(dtot-path.d[1])+hrs*path.d[1])/dtot;
    double hobs = HHi;
    double alpha_obt = HHi/path.d[1];
    double alpha_obr = HHi/(dtot-path.d[1]);

    //get max values
    for(int i=2; i<path.length-1;i++){
        HHi = path.h[i]-(hts*(dtot-path.d[i])+hrs*path.d[i])/dtot; //Eq 165d
        hobs = std::max(hobs, HHi);                                //Eq 165a
        alpha_obt = std::max(alpha_obt,HHi/path.d[i]);             //Eq 165b
        alpha_obr = std::max(alpha_obr,HHi/(dtot-path.d[i]));      //Eq 165c
    }

    double gt = alpha_obt/(alpha_obt+alpha_obr); //Eq 166e
    double gr = alpha_obr/(alpha_obt+alpha_obr); //Eq 166f

    double hstp = hst; //Eq 166a
    double hsrp = hsr; //Eq 166b
    if(hobs>0){
        hstp-=hobs*gt; //Eq 166c
        hsrp-=hobs*gr; //Eq 166d
    }
    
    SmoothEarthResults results = SmoothEarthResults();
    results.hst = hst;
    results.hsr = hsr;
    results.hstd = std::min(path.h.front(), hstp); //Eq 167 a,b
    results.hsrd = std::min(path.h.back(), hsrp); //Eq 167 c,d

    //horizon elevation angle and distance
    double theta_max;
    int kindex=0;
    double theta = 1000*std::atan((path.h[1]-hts)/1000*path.d[1])-path.d[1]/(2*ae);
    theta_max = theta;
    for(int i=2; i<path.length-1;i++){
        theta = 1000*std::atan((path.h[i]-hts)/1000*path.d[i])-path.d[i]/(2*ae);
        if(theta>theta_max){//in case of ties, favor the point closer to tx
            theta_max = theta;
            kindex = i;
        }
    }

    for(int i=2; i<path.length-1;i++){
        theta = 1000*std::atan((path.h[i]-hts)/1000*path.d[i])-path.d[i]/(2*ae);
        if(theta>=theta_max){//in case of ties, favor the point closer to rx
            theta_max = theta;
            kindex = i;
        }
    }

    return results;
}
struct SmoothEarthResults{
    double hst;
    double hsr;
    double hstd;
    double hsrd;
    double hte;
    double hre;
    double hm;
    double dlt;
    double dlr;
    double theta_t;
    double theta_r;
    double theta_tot;
    int pathtype;
};*/