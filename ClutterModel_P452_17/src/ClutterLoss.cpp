#include "ClutterModel_P452_17/ClutterLoss.h"
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <tuple>
#include <iostream>

namespace{
    const std::vector<ClutterModel::ClutterNominalHeight_m_andDistance_km> ClutterTable = {
        {0.0,0.0},
        {4.0,0.1},{4.0,0.1},{4.0,0.1},{4.0,0.1},{4.0,0.1},
        {5.0,0.07},
        {15.0,0.05},{15.0,0.05},{15.0,0.05},
        {20.0,0.05},{20.0,0.05},
        {20.0,0.03},
        {9.0,0.025},
        {12.0,0.02},
        {20.0,0.02},
        {25.0,0.02},
        {35.0,0.02},
        {20.0,0.05}
    };
}

ClutterModel::ClutterNominalHeight_m_andDistance_km ClutterModel::fetchNominalClutterValues(const ClutterType& clutterType){
    return ClutterTable.at(static_cast<int>(clutterType));
}

//WARNING ignoring site shielding for now
ClutterModel::ClutterResults ClutterModel::calculateClutterModel(const double& freq_GHz, const PathProfile::Path& path, 
        const double& height_tx_m, const double& height_rx_m, const ClutterType& tx_clutterType, 
        const ClutterType& rx_clutterType){

        
    const auto [tx_clutter_height_m,tx_clutter_dist_km] = fetchNominalClutterValues(tx_clutterType);
    const auto [rx_clutter_height_m,rx_clutter_dist_km] = fetchNominalClutterValues(rx_clutterType);

    uint16_t index1 = 0;
    uint16_t index2 = path.size();//sentinel index. the last valid value is right before this
    double hg_height_tx_m = height_tx_m;
    double hg_height_rx_m = height_rx_m;
    double tx_clutterLoss_dB = 0;
    double rx_clutterLoss_dB = 0;

    //make sure clutter model is applicable (clutter higher than antenna height)
    if(tx_clutter_height_m>height_tx_m){
        const double Ffc = 0.25+0.375*(1+std::tanh(7.5*(freq_GHz-0.5))); //Eq 57a
        tx_clutterLoss_dB = 10.25*Ffc*std::exp(-tx_clutter_dist_km)*(1-std::tanh(6*(height_tx_m/tx_clutter_height_m-0.625)))-0.33; //Eq 57

        //path length correction
        auto it = std::find_if(path.cbegin(),path.cend(),
            [tx_clutter_dist_km](PathProfile::ProfilePoint point){return point.d_km>=tx_clutter_dist_km;});
        if(it!=path.cend()){
            index1 = it-path.cbegin();
        }
        else{
            index1 = path.size();
        }
        hg_height_tx_m = tx_clutter_height_m;

    }

    if(rx_clutter_height_m>height_rx_m){
        const double Ffc = 0.25+0.375*(1+std::tanh(7.5*(freq_GHz-0.5))); //Eq 57a
        rx_clutterLoss_dB = 10.25*Ffc*std::exp(-rx_clutter_dist_km)*(1-std::tanh(6*(height_rx_m/rx_clutter_height_m-0.625)))-0.33; //Eq 57

        //path length correction
        const double rx_clutter_loc = path.back().d_km-rx_clutter_dist_km;
        auto it = std::find_if(path.cbegin(),path.cend(),
            [rx_clutter_loc](PathProfile::ProfilePoint point){return point.d_km>rx_clutter_loc;});
        if(it!=path.cend()){
            index2 = it-path.cbegin();
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
    const double offset = (*(path.begin()+index1)).d_km;
    //loop through middle segment of path
	ClutterResults resultObj;    
    PathProfile::ProfilePoint point;
    for(auto cit = path.cbegin()+index1; cit<path.cbegin()+index2; ++cit){
        point = *cit;
        resultObj.modifiedPath.push_back(PathProfile::ProfilePoint(point.d_km-offset, point.h_asl_m, point.zone));
    }
    resultObj.modifiedHeights_m = ITUR_P452::TxRxPair{hg_height_tx_m,hg_height_rx_m};
    resultObj.clutterLoss_dB = ITUR_P452::TxRxPair{tx_clutterLoss_dB,rx_clutterLoss_dB};

    return resultObj;
}