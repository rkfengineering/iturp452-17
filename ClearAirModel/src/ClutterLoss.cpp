#include "ClearAirModel/ClutterLoss.h"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <iostream>

ClearAirModel::ClutterLoss::ClutterLoss(const double& freq_GHz, const PathProfile::Path& path, 
        const double& height_tx_m, const double& height_rx_m, const ClutterType& tx_clutterType, 
        const ClutterType& rx_clutterType):
        m_freq_GHz{freq_GHz}, m_path{path}, m_height_tx_m{height_tx_m}, m_height_rx_m{height_rx_m}{

    populateNominalClutterValues(tx_clutterType, rx_clutterType);
    calcClutterLoss_dB();        
}

void ClearAirModel::ClutterLoss::populateNominalClutterValues(const ClutterType& tx_clutterType, const ClutterType& rx_clutterType){
    std::tie(m_tx_clutter_height_m,m_tx_clutter_dist_km) = ClutterTable.at(static_cast<int>(tx_clutterType));
    std::tie(m_rx_clutter_height_m,m_rx_clutter_dist_km) = ClutterTable.at(static_cast<int>(rx_clutterType));
}


//WARNING ignoring site shielding for now

void ClearAirModel::ClutterLoss::calcClutterLoss_dB(){
    u_int16_t index1 = 0;
    u_int16_t index2 = m_path.size();//sentinel index. the last valid value is right before this
    double hg_height_tx_m = m_height_tx_m;
    double hg_height_rx_m = m_height_rx_m;
    double tx_clutterLoss_dB = 0;
    double rx_clutterLoss_dB = 0;

    //make sure clutter model is applicable (clutter higher than antenna height)
    if(m_tx_clutter_height_m>m_height_tx_m){
        const double Ffc = 0.25+0.374*(1+std::tanh(7.5*(m_freq_GHz-0.5))); //Eq 57a
        tx_clutterLoss_dB = 10.25*Ffc*std::exp(-m_tx_clutter_dist_km)*(1-std::tanh(6*(m_height_tx_m/m_tx_clutter_height_m-0.625)))-0.33; //Eq 57

        //path length correction
        auto it = std::find_if(m_path.cbegin(),m_path.cend(),
            [this](PathProfile::ProfilePoint point){return point.d_km>=m_tx_clutter_dist_km;});
        if(it!=m_path.cend()){
            index1 = it-m_path.cbegin();
        }
        else{
            index1 = m_path.size();
        }
        hg_height_tx_m = m_tx_clutter_height_m;

    }

    if(m_rx_clutter_height_m>m_height_rx_m){
        const double Ffc = 0.25+0.374*(1+std::tanh(7.5*(m_freq_GHz-0.5))); //Eq 57a
        rx_clutterLoss_dB = 10.25*Ffc*std::exp(-m_rx_clutter_dist_km)*(1-std::tanh(6*(m_height_rx_m/m_rx_clutter_height_m-0.625)))-0.33; //Eq 57

        //path length correction
        const double rx_clutter_loc = m_path.back().d_km-m_rx_clutter_dist_km;
        auto it = std::find_if(m_path.cbegin(),m_path.cend(),
            [rx_clutter_loc](PathProfile::ProfilePoint point){return point.d_km>rx_clutter_loc;});
        if(it!=m_path.cend()){
            index2 = it-m_path.cbegin();
        }
        else{
            index2 = 0;
        }
        hg_height_rx_m = m_rx_clutter_height_m;
    }

    //error check
    //if (index2-index1<4){
        //std::cerr<<"sum of clutter nominal distances is larger than the path length"
    //}

    //TODO this can probably be optimized instead of just copying
    const double offset = (*(m_path.begin()+index1)).d_km;
    //loop through middle segment of path
    PathProfile::ProfilePoint point;
    for(auto cit = m_path.cbegin()+index1; cit<m_path.cbegin()+index2; ++cit){
        point = *cit;
        m_mod_path.push_back(PathProfile::ProfilePoint(point.d_km-offset, point.h_asl_m, point.zone));
    }
    m_hg_height_m = ClearAirModel::TxRxPair{hg_height_tx_m,hg_height_rx_m};
    m_clutterLoss_dB = ClearAirModel::TxRxPair{tx_clutterLoss_dB,rx_clutterLoss_dB};
}