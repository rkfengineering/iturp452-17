#include "ITUR_P452/PathProfile.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdint>

//constructors
PathProfile::ProfilePoint::ProfilePoint(){
}
PathProfile::ProfilePoint::ProfilePoint(double distance_km, double height_asl_m):
    d_km{distance_km}, h_asl_m{height_asl_m}{
}
PathProfile::ProfilePoint::ProfilePoint(double distance_km, double height_asl_m, ZoneType zonetype):
    d_km{distance_km}, h_asl_m{height_asl_m}, zone{zonetype}{
}

//Default Constructor
PathProfile::Path::Path(){
}

//TODO we can add an assumption that height==0 implies over sea. otherwise land. 
//TODO look for a more standard approach to interfacing with csvs
PathProfile::Path::Path(std::string csvPath){
    std::ifstream file;
    file.open(csvPath);
    std::string line;
    std::getline(file,line);//throw away first header line

    std::stringstream ss(line);
    std::string val;
    uint16_t n = 0;
    while(std::getline(ss,val,',')){
        n+=1;
    }

    while(std::getline(file,line)){
        std::stringstream ss(line);
        std::string val;
        PathProfile::ProfilePoint point;

        for(uint16_t i=0; i<n; i++){
            std::getline(ss,val,',');
            //row 1 is the distance value
            if(i==0){
                point.d_km = std::stod(val);
            }
            //row 2 is the height value (m)
            else if(i==1){
                point.h_asl_m = std::stod(val);
            }
            //Row 3 is a label of the zone type, can be skipped
            //Row 4 is an integer representing zone type
            else if(i==3){
                point.zone = static_cast<PathProfile::ZoneType>(std::stoi(val));
            }
        }
        push_back(point);
    }
}

double PathProfile::Path::calcFracOverSea() const{
    double sea_dist = 0;
    PathProfile::ProfilePoint lastPoint = *cbegin();
    for(auto cit = cbegin()+1; cit<cend(); ++cit){
        //add full interval (both points are sea type)
        if(cit->zone==PathProfile::ZoneType::Sea && lastPoint.zone==PathProfile::ZoneType::Sea){
            sea_dist += cit->d_km - lastPoint.d_km;
        }
        //add half interval (transition from sea to land or land to sea)
        else if (cit->zone==PathProfile::ZoneType::Sea || lastPoint.zone==PathProfile::ZoneType::Sea){
            sea_dist += (cit->d_km - lastPoint.d_km)/2.0;
        }
        lastPoint = *cit;
    }
    const double full_dist = back().d_km;
    return sea_dist/full_dist;
}

//WARNING this function only works if the zone types are populated correctly (not checked)
double PathProfile::Path::calcTimePercentBeta0(const double& centerLatitude_deg) const{
    //get longest contiguous land and inland segments 

    double longestLand=0;
    double longestInland=0;
    double currentLand=0;
    double currentInland=0;
    
    PathProfile::ProfilePoint lastPoint = *cbegin();

    for(auto cit = cbegin()+1; cit<cend(); ++cit){
        //add full interval (both points are land type)
        if(cit->zone!=PathProfile::ZoneType::Sea && lastPoint.zone!=PathProfile::ZoneType::Sea){
            currentLand += cit->d_km - lastPoint.d_km;
        }
        //or add half interval (transition from sea to land or land to sea)
        else if (cit->zone==PathProfile::ZoneType::Sea || lastPoint.zone==PathProfile::ZoneType::Sea){
            currentLand += (cit->d_km - lastPoint.d_km)/2.0;
            //reset 
            if (cit->zone==PathProfile::ZoneType::Sea){
                longestLand = std::max(longestLand,currentLand);
                currentLand = 0;
            }
        }

        //add full interval (both points are inland type)
        if(cit->zone==PathProfile::ZoneType::Inland && lastPoint.zone==PathProfile::ZoneType::Inland){
            currentInland += cit->d_km - lastPoint.d_km;
        }
        //or add half interval (transition to or from inland)
        else if (cit->zone==PathProfile::ZoneType::Inland || lastPoint.zone==PathProfile::ZoneType::Inland){
            currentInland += (cit->d_km - lastPoint.d_km)/2.0;
            //reset 
            if (cit->zone!=PathProfile::ZoneType::Inland){
                longestInland = std::max(longestInland,currentInland);
                currentInland = 0;
            }
        }
        lastPoint = *cit;
    }
    //max values only update on transition out of zone. 
    //check if the longest contiguous zone is at the end of the path
    longestLand = std::max(longestLand,currentLand);
    longestInland = std::max(longestInland,currentInland);

    //calculate beta0

    //Equation 3a
    const double tau = 1.0-std::exp(-(4.12*1e-4*std::pow(longestInland,2.41)));
    //Equation 3
    const double mu1a = std::pow(10.0,-longestLand/(16.0-6.6*tau));
    const double mu1b = std::pow(10.0,-5*(0.496+0.354*tau));
    //limit to mu<=1
    const double mu1 = std::min(std::pow(mu1a+mu1b, 0.2),1.0);

    const double abs_phi = std::abs(centerLatitude_deg);
    if (abs_phi<=70){
        //Equation 4
        const double mu4 = std::pow(10.0,(-0.935+0.0176*abs_phi)*std::log10(mu1));
        //Equation 2
        return std::pow(10.0, -0.015*abs_phi + 1.67)*mu1*mu4;
    }
    else{
        //Equation 4
        const double mu4 = std::pow(10.0,0.3*std::log10(mu1));
        //Equation 2
        return 4.17*mu1*mu4;
    }
}

//TODO refactor code. calc beta0 needs inland and non-sea
// ducting model needs inland only
double PathProfile::Path::calcLongestContiguousInlandDistance_km() const{
    //get longest contiguous land and inland segments 

    double longestInland=0;
    double currentInland=0;
    
    PathProfile::ProfilePoint lastPoint = *cbegin();

    for(auto cit = cbegin()+1; cit<cend(); ++cit){

        //add full interval (both points are inland type)
        if(cit->zone==PathProfile::ZoneType::Inland && lastPoint.zone==PathProfile::ZoneType::Inland){
            currentInland += cit->d_km - lastPoint.d_km;
        }
        //or add half interval (transition to or from inland)
        else if (cit->zone==PathProfile::ZoneType::Inland || lastPoint.zone==PathProfile::ZoneType::Inland){
            currentInland += (cit->d_km - lastPoint.d_km)/2.0;
            //reset 
            if (cit->zone!=PathProfile::ZoneType::Inland){
                longestInland = std::max(longestInland,currentInland);
                currentInland = 0;
            }
        }
        lastPoint = *cit;
    }
    //max values only update on transition out of zone. 
    //check if the longest contiguous zone is at the end of the path
    longestInland = std::max(longestInland,currentInland);
    return longestInland;
}