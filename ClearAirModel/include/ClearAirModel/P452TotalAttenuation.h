#ifndef P452_TOTAL_ATTENUATION_H
#define P452_TOTAL_ATTENUATION_H

#include "PathProfile.h"
#include "ClearAirModelHelpers.h"
#include "BasicProp.h"
#include "DiffractionLoss.h"
#include "TropoScatter.h"
#include "AnomolousProp.h"
#include "ClutterLoss.h"
#include "Common/GeodeticCoord.h"
#include "Common/Enumerations.h"

namespace ClearAirModel{

//Section 4.6  Basic transmission loss between the two stations
class p452_TotalAttenuation {
public:

    /// @brief Calculates Basic transmission loss (dB), not exceeded for the required annual percentage time, p,
    /// @param freq_GHz             Frequency (GHz)
    /// @param p_percent            Percentage of time not exceeded (%), 0<p<=50
    /// @param path                 Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param height_tx_m          Tx Antenna height (m)
    /// @param height_rx_m          Rx Antenna height (m)
    /// @param centerLatitude_deg   The latitude of the path center point (deg) 
    /// @param txHorizonGain_dBi    Tx Antenna directional gain towards the horizon along the path (dB)
    /// @param rxHorizonGain_dBi    Rx Antenna directional gain towards the horizon along the path (dB)
    /// @param pol                  Polarization type (horizontal or vertical)
    /// @param dist_coast_tx_km 
    /// @param dist_coast_rx_km 
    /// @param deltaN 
    /// @param surfaceRefractivity 
    /// @param temp_K 
    /// @param dryPressure_hPa 
    /// @param tx_clutterType
    /// @param rx_clutterType
    /// @return 
    p452_TotalAttenuation(const double& freq_GHz, const double& p_percent, const PathProfile::Path path, 
            const double& height_tx_m, const double& height_rx_m, const double& centerLatitude_deg, const double& txHorizonGain_dBi, 
            const double& rxHorizonGain_dBi, const Enumerations::PolarizationType& pol, const double& dist_coast_tx_km, 
            const double& dist_coast_rx_km, const double& deltaN, const double& surfaceRefractivity,
            const double& temp_K, const double& dryPressure_hPa, const ClutterType& tx_clutterType, 
            const ClutterType& rx_clutterType);

    inline double getTotalTransmissionLoss_dB() const{ return m_totalTransmissionLoss_dB;}

private:
    //common direct inputs
    const double& m_freq_GHz;
    const double& m_p_percent;

    //height gain model variables
    PathProfile::Path m_mod_path;
    double m_height_tx_asl_m; 
    double m_height_rx_asl_m;
    double m_d_tot_km;

    //intermediate inputs
    ClearAirModel::HorizonAnglesAndDistances m_HorizonVals;
    double m_b0_percent;
    double m_fracOverSea;
    double m_effEarthRadius_med_km;
    double m_horizonAnglesAndDistances;

    //submodel outputs
    double m_freeSpaceWithGasLoss_dB;
    double m_basicTransmissionLoss_p_percent_dB;
    double m_basicTransmissionLoss_b0_percent_dB;
    double m_diffractionLoss_p_percent_dB;
    double m_diffractionLoss_median_dB;
    double m_anomolousPropagationLoss_dB;
    double m_tropoScatterLoss_dB;
    double m_tx_clutterLoss_dB;
    double m_rx_clutterLoss_dB;

    //final value
    double m_totalTransmissionLoss_dB;

    //Apply clutter/height gain model and calculate path parameters
    void populatePathParameters(const PathProfile::Path& path, const double& deltaN, const double& centerLatitude_deg,
            const double& height_tx_m, const double& height_rx_m, const ClutterType& tx_clutterType, 
            const ClutterType& rx_clutterType);

    //TODO replace DN with median effective earth radius as input for diffraction model
    void calculateSubModels(const double& temp_K, const double& dryPressure_hPa, const double deltaN,
        const double& dist_coast_tx_km, const double& dist_coast_rx_km, const double& seaLevelSurfaceRefractivity, 
        const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi, const Enumerations::PolarizationType& pol);

    double calcP452TotalAttenuation();

    double calcSlopeInterpolationParameter(const PathProfile::Path path, const double effEarthRadius_med_km,
            const double& height_tx_asl_m,const double& height_rx_asl_m);
            
    double calcPathBlendingInterpolationParameter(const double& d_tot_km);
};//end class p452_TotalAttenuation

} //end namespace ClearAirModel

#endif /* P452_TOTAL_ATTENUATION_H */