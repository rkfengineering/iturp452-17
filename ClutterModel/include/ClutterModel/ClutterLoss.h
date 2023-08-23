#ifndef CLUTTER_LOSS_H
#define CLUTTER_LOSS_H

#include "MainModel/PathProfile.h"
#include "MainModel/Helpers.h"
#include <vector>

namespace ClutterModel {
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

struct ClutterResults{
    //distances (km), heights (asl)(m), and zone types of the profile points in the height gain model
    PathProfile::Path modifiedPath;
    //Tx,Rx Antenna height above ground level (m) in the height gain model
    ITUR_P452::TxRxPair modifiedHeights_m;
    //Additional clutter shielding losses at Tx,Rx (dB)
    ITUR_P452::TxRxPair clutterLoss_dB;
};

/// @brief Main entry point to execute height gain model calculations 
///        (returns modified path, modified antenna heights, additional clutter losses)
/// @param freq_GHz             Transmitting Frequency (GHz) 
/// @param path                 distances (km), heights (asl_m), and zone types of the profile points
/// @param height_tx_m          Tx Antenna center height above ground level (m)
/// @param height_rx_m          Rx Antenna center height above ground level (m)
/// @param tx_clutterType       Clutter Category Type at tx 
/// @param rx_clutterType       Clutter Category Type at rx 
ClutterResults calculateClutterModel(const double& freq_GHz, const PathProfile::Path& path, 
        const double& height_tx_m, const double& height_rx_m, const ClutterType& tx_clutterType, 
        const ClutterType& rx_clutterType);
    
/// @brief fetch Nominal Height (m) and Distance(km) values for clutter type using table 4
/// @param tx_clutterType Clutter Type at Tx
/// @param rx_clutterType Clutter Type at Rx
ClutterNominalHeight_m_andDistance_km fetchNominalClutterValues(const ClutterType& clutterType);

};//end namespace ClutterModel
#endif /* CLUTTER_LOSS_H */