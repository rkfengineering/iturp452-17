#include "ClutterLoss.h"
#include <cmath>
#include <algorithm>

//WARNING ignoring site shielding for now

ClutterLoss::ClutterLossResults ClutterLoss::clutterLoss_corr(const double& freq_GHz, 
            const PathProfile::Path& path, const double& height_tx_m, const double& height_rx_m, 
            const double& tx_clutter_height_m, const double& rx_clutter_height_m, const double& tx_clutter_dist_km,
            const double& rx_clutter_dist_km){

    int index1 = 0;
    int index2 = path.size()-1;
    double hg_height_tx_m = height_tx_m;
    double hg_height_rx_m = height_rx_m;
    double tx_clutterLoss_dB = 0;
    double rx_clutterLoss_dB = 0;

    //make sure clutter model is applicable (clutter higher than antenna height)
    if(tx_clutter_height_m>height_tx_m){
        const double Ffc = 0.25+0.374*(1+std::tanh(7.5*(freq_GHz-0.5))); //Eq 57a
        tx_clutterLoss_dB = 10.25*Ffc*std::exp(-tx_clutter_dist_km)*(1-std::tanh(6*(height_tx_m/tx_clutter_height_m-0.625)))-0.33; //Eq 57

        //path length correction
        auto it = std::find_if(path.begin(),path.end(),
            [tx_clutter_dist_km](PathProfile::ProfilePoint point){return point.d_km>=tx_clutter_dist_km;});
        if(it!=path.end()){
            index1 = it-path.begin();
        }
        else{
            index1 = path.size()-1;
        }
        hg_height_tx_m = tx_clutter_height_m;

    }

    if(rx_clutter_height_m>height_rx_m){
        const double Ffc = 0.25+0.374*(1+std::tanh(7.5*(freq_GHz-0.5))); //Eq 57a
        rx_clutterLoss_dB = 10.25*Ffc*std::exp(-rx_clutter_dist_km)*(1-std::tanh(6*(height_rx_m/rx_clutter_height_m-0.625)))-0.33; //Eq 57

        //path length correction
        const double rx_clutter_loc = path.back().d_km-rx_clutter_dist_km;
        auto it = std::find_if(path.begin(),path.end(),
            [rx_clutter_loc](PathProfile::ProfilePoint point){return point.d_km>rx_clutter_loc;});
        if(it!=path.end()){
            index2 = it-path.begin();
        }
        else{
            index2 = 0;
        }
        hg_height_rx_m = rx_clutter_height_m;
    }

    //error check
    //if (index2-index1<4){
        //std::cerr<<"sum of clutter nominal distances is larger than the path length"
    //}

    //TODO this can probably be optimized instead of just copying
    PathProfile::Path modified_path;
    const double offset = (*(path.begin()+index1)).d_km;
    //loop through middle segment of path
    PathProfile::ProfilePoint point;
    for(auto cit = path.cbegin()+index1; cit<path.cbegin()+index2; ++cit){
        point = *cit;
        modified_path.push_back(PathProfile::ProfilePoint(point.d_km-offset, point.h_masl, point.zone));
    }

    ClutterLossResults results;
    results.path = modified_path;
    results.hg_height_tx_m = hg_height_tx_m;
    results.hg_height_rx_m = hg_height_rx_m;
    results.tx_clutterLoss_dB = tx_clutterLoss_dB;
    results.rx_clutterLoss_dB = rx_clutterLoss_dB;

    return results;
}