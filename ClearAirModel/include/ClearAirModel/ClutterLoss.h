#ifndef CLUTTER_LOSS_H
#define CLUTTER_LOSS_H

#include "PathProfile.h"
#include "ClearAirModelHelpers.h"
#include <vector>

namespace ClearAirModel{

//Section 4.5 Additional clutter losses

//Pair for tx (first) and rx (second) values
using ClutterNominalHeight_m_andDistance_km = std::pair<double,double>;
/// Specifies Clutter Category Type
enum ClutterType {
    NoClutter = 0,
    HighCropFields,
    ParkLand,
    IrregularlySpacedSparseTrees,
    Orchard_RegularlySpaced,
    SparseHouses,
    VillageCentre,
    DeciduousTrees_IrregularlySpaced,
    DeciduousTrees_RegularlySpaced,
    MixedTreeForest,
    ConiferousTrees_IrregularlySpaced,
    ConiferousTrees_RegularlySpaced,
    TropicalRainForest,
    Suburban,
    DenseSuburban,
    Urban,
    DenseUrban,
    HighRiseUrban,
    IndustrialZone
};

//TODO find a better place to put this information
static const std::vector<ClutterNominalHeight_m_andDistance_km> ClutterTable = {
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

class ClutterLoss {
public:
    
    /// @brief Compute the height-gain correlation from Section 4.5.4. Uses TABLE 4 for Nominal height/distance
    /// @param freq_GHz             Transmitting Frequency (GHz) 
    /// @param path                 distances (km), heights (asl_m), and zone types of the profile points
    /// @param height_tx_m          Tx Antenna center height above ground level (m)
    /// @param height_rx_m          Rx Antenna center height above ground level (m)
    /// @param tx_clutterType       Clutter Category at tx 
    /// @param rx_clutterType       Clutter Category at rx 
    ClutterLoss(const double& freq_GHz, const PathProfile::Path& path, 
            const double& height_tx_m, const double& height_rx_m, const ClutterType& tx_clutterType, 
            const ClutterType& rx_clutterType);

    PathProfile::Path getModifiedPath() const {return m_mod_path;}
    ClearAirModel::TxRxPair getHeightGainModelHeights_m() const{ return m_hg_height_m; }
    ClearAirModel::TxRxPair getClutterLoss_dB() const{ return m_clutterLoss_dB; }

private:
    const double& m_freq_GHz;
    const PathProfile::Path& m_path;
    const double& m_height_tx_m;       //Tx Antenna center height above ground level (m)
    const double& m_height_rx_m;       //Rx Antenna center height above ground level (m)    

    double m_tx_clutter_height_m; //Values from lookup table
    double m_rx_clutter_height_m; 
    double m_tx_clutter_dist_km;
    double m_rx_clutter_dist_km;

    /// distances (km), heights (asl_m), and zone types of the profile points in the height gain model
    PathProfile::Path m_mod_path; 
    /// Tx,Rx Antenna center height above ground level (m) in the height gain model
    ClearAirModel::TxRxPair m_hg_height_m;
    /// Additional clutter shielding losses at Tx,Rx (dB)
    ClearAirModel::TxRxPair m_clutterLoss_dB;
    
    /// @brief Populates Nominal Height (m) and Distance(km) values for clutter types using table
    /// @param tx_clutterType Clutter Type at Tx
    /// @param rx_clutterType Clutter Type at Rx
    void populateNominalClutterValues(const ClutterType& tx_clutterType, const ClutterType& rx_clutterType);

    /// @brief Calculates clutter losses, height gain model antenna heights, and height gain model modified path
    void calcClutterLoss_dB();

};//end class ClutterLoss

} //end namespace ClearAirModel
#endif /* CLUTTER_LOSS_H */