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

