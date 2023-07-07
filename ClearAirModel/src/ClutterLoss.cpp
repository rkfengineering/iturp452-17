#include "ClutterLoss.h"
#include <cmath>
#include <algorithm>

//WARNING ignoring site shielding for now

ClutterLoss::ClutterLossResults ClutterLoss::closs_corr(const double freq_GHz, 
            const ProfilePath& path, const double htg, const double hrg, 
            const double ha_t, const double ha_r, const double dk_t, const double dk_r){

    int index1 = 0;
    int index2 = path.length-1;
    double htgc = htg;
    double hrgc = hrg;
    double Aht = 0;
    double Ahr = 0;
    double ha = ha_t;
    double dk = dk_t;

    //make sure clutter model is applicable (clutter higher than antenna height)
    if(ha>htg){
        double Ffc = 0.25+0.374*(1+std::tanh(7.5*(freq_GHz-0.5))); //Eq 57a
        Aht = 10.25*Ffc*std::exp(-dk)*(1-std::tanh(6*(htg/ha-0.625)))-0.33; //Eq 57

        //path length correction
        auto it = std::find_if(path.d.begin(),path.d.end(),[dk](double point){return point>=dk;});
        if(it!=path.d.end()){
            index1 = it-path.d.begin();
        }
        else{
            index1 = path.length-1;
        }
        htgc = ha_t;

    }

    ha = ha_r;
    dk = dk_r;

    if(ha>hrg){
        double Ffc = 0.25+0.374*(1+std::tanh(7.5*(freq_GHz-0.5))); //Eq 57a
        Ahr = 10.25*Ffc*std::exp(-dk)*(1-std::tanh(6*(hrg/ha-0.625)))-0.33; //Eq 57

        //path length correction
        double rx_clutter_loc = path.d.back()-dk;
        auto it = std::find_if(path.d.begin(),path.d.end(),[rx_clutter_loc](double point){return point>rx_clutter_loc;});
        if(it!=path.d.end()){
            index2 = it-path.d.begin();
        }
        else{
            index2 = 0;
        }
        hrgc = ha_r;
    }

    //error check
    //if (index2-index1<4){
        //std::cerr<<"sum of clutter nominal distances is larger than the path length"
    //}

    //TODO this can probably be optimized instead of just copying
    ProfilePath modified_path;
    modified_path.d = std::vector<double>(path.d.begin()+index1, path.d.begin()+index2);
    double offset = *(path.d.begin()+index1);
    std::for_each(modified_path.d.begin(), modified_path.d.end(), [offset](double& x){x -= offset;});

    modified_path.h = std::vector<double>(path.h.begin()+index1, path.h.begin()+index2);
    modified_path.zone = std::vector<int>(path.zone.begin()+index1, path.zone.begin()+index2);
    modified_path.length = index2-index1;

    ClutterLossResults results;
    results.path = modified_path;
    results.htgc = htgc;
    results.hrgc = hrgc;
    results.Aht = Aht;
    results.Ahr = Ahr;

    return results;
}