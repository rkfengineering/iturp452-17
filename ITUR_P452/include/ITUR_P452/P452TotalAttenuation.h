#ifndef P452_TOTAL_ATTENUATION_H
#define P452_TOTAL_ATTENUATION_H

#include "PathProfile.h"
#include "Helpers.h"
#include "BasicProp.h"
#include "DiffractionLoss.h"
#include "TropoScatter.h"
#include "AnomalousProp.h"
#include "ClutterModel_P452_17/ClutterLoss.h"
#include "Common/GeodeticCoord.h"
#include "Common/Enumerations.h"

namespace ITUR_P452{

//Section 4.6  Basic transmission loss between the two stations
class TotalClearAirAttenuation {
public:
    /// @brief Calculates Basic transmission loss (dB), not exceeded for the required annual percentage time, p,
    /// @param freq_GHz             Frequency (GHz)
    /// @param p_percent            Percentage of time not exceeded (%), 0<p<=50
    /// @param path_TxToRx          Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param height_tx_m          Tx Antenna height (m)
    /// @param height_rx_m          Rx Antenna height (m)
    /// @param centerLatitude_deg   The latitude of the path center point (deg) 
    /// @param txHorizonGain_dBi    Tx Antenna directional gain towards the horizon along the path (dB)
    /// @param rxHorizonGain_dBi    Rx Antenna directional gain towards the horizon along the path (dB)
    /// @param pol                  Polarization type (horizontal or vertical)
    /// @param dist_coast_tx_km     Distance over land from Tx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param dist_coast_rx_km     Distance over land from Rx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param deltaN               Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (N-Units/km) 
    /// @param surfaceRefractivity  Sea Level Surface Refractivity (N0) (N-Units)
    /// @param temp_K               Temperature (K)
    /// @param dryPressure_hPa      Dry air pressure (hPa)
    /// @param tx_clutterType       Clutter Category Type at Tx 
    /// @param rx_clutterType       Clutter Category Type at Rx 
    TotalClearAirAttenuation(const double& freq_GHz, const double& p_percent, const PathProfile::Path path_TxToRx, 
            const double& height_tx_m, const double& height_rx_m, const double& centerLatitude_deg, const double& txHorizonGain_dBi, 
            const double& rxHorizonGain_dBi, const Enumerations::PolarizationType& pol, const double& dist_coast_tx_km, 
            const double& dist_coast_rx_km, const double& deltaN, const double& surfaceRefractivity,
            const double& temp_K, const double& dryPressure_hPa, const ClutterModel::ClutterType& tx_clutterType, 
            const ClutterModel::ClutterType& rx_clutterType);

    /// @brief Combine submodel results according to the method in Section 4.6
    /// @return total transmission loss for clear air conditions
    double calcTotalClearAirAttenuation() const;

private:
    //common direct inputs
    const double& m_freq_GHz;       //Frequency (GHz)
    const double& m_p_percent;      //Percentage of time not exceeded (%), 0<p<=50

    //height gain model variables
    PathProfile::Path m_mod_path;   //distances (km), heights (asl)(m), and zone types of the profile points in the height gain model
    double m_height_tx_asl_m;       //Tx Antenna center height above ground level (m)
    double m_height_rx_asl_m;       //Rx Antenna center height above ground level (m)
    double m_d_tot_km;              //Great Circle Distance between Tx and Rx antennas along modified path (km)

    //intermediate inputs

    //Tx and Rx Horizon Elevation Angles (mrad) and Tx and Rx Horizon Distances (km)
    ITUR_P452::HorizonAnglesAndDistances m_HorizonVals;
    double m_b0_percent;                   //Time percentage that the refractivity gradient exceeds 100 N-Units/km
    double m_fracOverSea;                  //Fraction of the path over sea
    double m_effEarthRadius_med_km;        //Median effective Earth's radius (km)

    //intermediate submodel outputs
    double m_freeSpaceWithGasLoss_dB;             //free space transmission loss with gas attenuation
    double m_basicTransmissionLoss_p_percent_dB;  //free space loss with gas atten and multipath focusing correction for p percent of time
    double m_basicTransmissionLoss_b0_percent_dB; //free space loss with gas atten and multipath focusing correction for b0 percent of time
    double m_diffractionLoss_p_percent_dB;        //Get diffraction Loss not exceeded for p percent of time (annual) (dB)
    double m_diffractionLoss_median_dB;           //Get diffraction Loss not exceeded for p=50% (dB)
    double m_anomalousPropagationLoss_dB;         //Transmission Loss with ducting and layer reflection (dB)
    double m_tropoScatterLoss_dB;                 //Loss due to troposcatter (dB)
    double m_tx_clutterLoss_dB;                   //loss associated with clutter shielding at tx
    double m_rx_clutterLoss_dB;                   //loss associated with clutter shielding at rx
    
    /// @brief Apply clutter/height gain model and calculate path parameters
    /// @param path                 Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param deltaN               Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (positive value) 
    /// @param centerLatitude_deg   The latitude of the path center point (deg) 
    /// @param height_tx_m          Tx Antenna height (m)
    /// @param height_rx_m          Rx Antenna height (m)
    /// @param tx_clutterType       Clutter Category Type at Tx
    /// @param rx_clutterType       Clutter Category Type at Rx
    void pre_calcPathParameters(const PathProfile::Path& path, const double& deltaN, const double& centerLatitude_deg,
            const double& height_tx_m, const double& height_rx_m, const ClutterModel::ClutterType& tx_clutterType, 
            const ClutterModel::ClutterType& rx_clutterType);

    //TODO replace DN with median effective earth radius as input for diffraction model

    /// @brief Create submodel objects and calculate submodel transmission loss outputs
    /// @param temp_K           Temperature (K)
    /// @param dryPressure_hPa  Dry air pressure (hPa)
    /// @param deltaN           Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (positive value) 
    /// @param dist_coast_tx_km Distance over land from Tx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param dist_coast_rx_km Distance over land from Rx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param seaLevelSurfaceRefractivity Sea Level Surface Refractivity (N0) (N-Units)
    /// @param txHorizonGain_dBi Tx Antenna directional gain towards the horizon along the path (dB)
    /// @param rxHorizonGain_dBi Rx Antenna directional gain towards the horizon along the path (dB)
    /// @param pol               Polarization type (horizontal or vertical)
    void calculateSubModels(const double& temp_K, const double& dryPressure_hPa, const double deltaN,
        const double& dist_coast_tx_km, const double& dist_coast_rx_km, const double& seaLevelSurfaceRefractivity, 
        const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi, const Enumerations::PolarizationType& pol);

    /// @brief calculate slope interpolation parameter used in Section 4.6
    /// @param path_TxToRx Profile path of distance (km) and height (asl) (m) points
    /// @param effEarthRadius_med_km Median effective earth radius (km)
    /// @param height_tx_asl_m Tx Antenna height above sea level (m)
    /// @param height_rx_asl_m Rx Antenna height above sea level (m)
    /// @return Slope Interpolation Parameter
    static double calcSlopeInterpolationParameter(const PathProfile::Path& path_TxToRx, const double& effEarthRadius_med_km,
            const double& height_tx_asl_m,const double& height_rx_asl_m);
    
    /// @brief calculate Path Blending interpolation parameter used in Section 4.6
    /// @param d_tot_km Distance between Tx and Rx antennas along great circle path (km)
    /// @return Path Blending Interpolation Parameter
    static double calcPathBlendingInterpolationParameter(const double& d_tot_km);
};//end class TotalClearAirAttenuation

} //end namespace ITUR_P452

#endif /* P452_TOTAL_ATTENUATION_H */