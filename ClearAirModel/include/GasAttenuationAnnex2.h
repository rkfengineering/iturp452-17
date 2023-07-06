#ifndef GAS_ATTENUATION_ANNEX2_H
#define GAS_ATTENUATION_ANNEX2_H

#include "DataLoader.h"
#include "Enumerations.h"
#include "GeodeticCoord.h"

class GasAttenuationAnnex2 {
public:
    GasAttenuationAnnex2() {}

    /// @brief Section 2.2 "Slant paths" of ITU-R P.676-12
    /// @param txLocation Transmitting location (lat/lon)
    /// @param freq_GHz Transmitting frequency (GHz)
    /// @param rho0_gm3 Surface water vapor density at the given location (gm^3)
    /// @param season Season in which the system is operating (winter/summer)
    /// @return Specific gaseous attenuation at a given transmitting location (dB/km)
    double calculateZenithAttenuation_dB(const GeodeticCoord& txLocation, const double& freq_GHz, 
                const double& rho0_gm3, const Enumerations::Season& season);

    double calculateEarthSpacePathAttenuation_dB(const double& height_km, const double& freq_GHz, 
                const double& elevationAngle_deg, const double& temp_K, const double& totalPressure_hPa, 
                const double& waterPressure_hPa, const double& integratedWaterVapor_kgm2);

    /// @brief 
    /// @param txLocation Transmitting location (lat/lon)
    /// @param freq_GHz Transmitting frequency (GHz)
    /// @param rho0_gm3 Surface water vapor density at the given location (g/m^3)
    /// @param integratedWaterVapor_kgm2 Integrated water vapor at the given location (kg/m^2)
    /// @param elevationAngle_deg Transmitting antenna's elevation angle (deg)
    /// @param season Season in which the system is operating (winter/summer)
    /// @return Attenuation along terrestrial path between transmitter and receiver due to gases (dB)
    double calculateEarthSpacePathAttenuationSeasonal_dB(const GeodeticCoord& txLocation, 
                const double& freq_GHz, const double& rho0_gm3, const double& integratedWaterVapor_kgm2, 
                const double& elevationAngle_deg, const Enumerations::Season& season);

    /// <summary>
    /// Section 2.2.1.2 "Inclined paths"
    /// </summary>
    /// <param name="txLocation">Transmitting location (lat/lon)</param>
    /// <param name="freq_GHz">Transmitting frequency (GHz)</param>
    /// <param name="avail">Desired availability (0-1)</param>
    /// <param name="elevationAngle_deg">Apparent elevation angle from the transmitting location to the space station (deg)</param>
    /// <returns>Attenuation along the slant path connecting the transmitting Earth station to the receiving space station (dB)</returns>
    double calculateEarthInclinedPathAttenuation_dB(const GeodeticCoord& txLocation, 
                const GeodeticCoord& rxLocation, const double& freq_GHz, const double& rho0_gm3, 
                const double& elevationAngle_deg, const Enumerations::Season& season);
};

#endif // GAS_ATTENUATION_ANNEX2_H