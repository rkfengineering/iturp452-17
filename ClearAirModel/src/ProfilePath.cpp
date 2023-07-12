#include "ProfilePath.h"
#include <fstream>
#include <sstream>

ProfilePath::ProfilePath(){
    d = std::vector<double>{};
    h = std::vector<double>{};
    zone = std::vector<int>{};
    length = 0;
}

//TODO look for a more standard approach to interfacing with csvs
ProfilePath::ProfilePath(std::string csvPath){
    std::ifstream file;
    file.open(csvPath);
    std::string line;
    std::getline(file,line);//throw away first header line

    std::stringstream ss(line);
    std::string val;
    int n = 0;
    while(std::getline(ss,val,',')){
        n+=1;
    }

    while(std::getline(file,line)){
        std::stringstream ss(line);
        std::string val;
        for(int i=0; i<n; i++){
            std::getline(ss,val,',');
            if(i==0){
                d.push_back(std::stod(val));
            }
            else if(i==1){
                h.push_back(std::stod(val));
            }
            else if(i==2){
                zone.push_back(std::stoi(val));
            }
        }
    }
    length = d.size();
}