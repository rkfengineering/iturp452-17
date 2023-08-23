#ifndef P452_H
#define P452_H

#include "MainModel/PathProfile.h"
#include "ClutterModel/ClutterLoss.h"
#include <vector>

//WARNING ITU_R P452 is not recommended for frequencies below 100 MHz (VHF band)
namespace P452 {

    /////////////////////////////
    // Main P452 Functions
    
    /// @brief Calculate total path loss for clear air conditions using ITU-R P.452-17 model, assuming summer season
    /// @param txHeight_m           Tx Antenna Height above terrain (m)
    /// @param rxHeight_m           Rx Antenna Height above terrain (m)
	/// @param elevationList_m      raw elevation list (meters above sea level) from tx to rx, total distance recommended <10,000 km
    /// @param stepDistance_km      distance between points in elevation list (km)
    /// @param midpoint_lat_deg     Latitude of midpoint in great circle path between tx and rx (deg)
    /// @param midpoint_lon_deg     Longitude of midpoint in great circle path between tx and rx (deg)
    /// @param freq_GHz             Frequency (GHz) Recommended between 0.1 GHz and 50 GHz
    /// @param timePercent          Required time percentage for which the calculated loss is not exceeded, 0<p<=50
    /// @param polariz              0 for Horizonatal Polarization, 1 for Vertical Polarization
    /// @param txHorizonGain_dBi    Tx Antenna directional gain towards the horizon along the path (dB)
    /// @param rxHorizonGain_dBi    Rx Antenna directional gain towards the horizon along the path (dB)
    /// @param txClutterType        Clutter Category Type at Tx 
    /// @param rxClutterType        Clutter Category Type at Rx 
	/// @return Path Loss (dB)
    double calculateP452Loss_dB(const double& txHeight_m, const double& rxHeight_m, 
            const std::vector<double>& elevationList_m, const double& stepDistance_km, 
            const double& midpoint_lat_deg, const double& midpoint_lon_deg,
            const double& freq_GHz, const double& timePercent, const int& polariz=0,
            const double& txHorizonGain_dBi=0, const double& rxHorizonGain_dBi=0,
            const ClutterModel::ClutterType& txClutterType=ClutterModel::ClutterType::NoClutter,
            const ClutterModel::ClutterType& rxClutterType=ClutterModel::ClutterType::NoClutter);

    /////////////////////////////
    // P452 Helper Functions

    /// @brief create path for ITU-R P.452-17 model assuming elevation==0 implies Sea zone type
	/// @param elevationList_m      raw elevation list (meters above sea level)
    /// @param stepDistance_km      distance between points in elevation list (km)
    /// @param out_path             return P452 path
	/// @param out_dist_coast_tx_km return distance from tx to the coast (km), 0 if at sea
	/// @param out_dist_coast_rx_km return distance from rx to the coast (km), 0 if at sea
    void createP452Path(const std::vector<double>& elevationList_m, const double& stepDistance_km,
        PathProfile::Path& out_path, double& out_dist_coast_tx_km, double& out_dist_coast_rx_km);

    /// Potentially Useful Functions:
    /// @brief create raw elevation list for ITU-R P.452-17 model using GdalRasterProcessor
    /// @brief create path for ITU-R P.452-17 model using land border data to classify zones through GdalVectorProcessor

} // end namespace P452
#endif /* P452_H */