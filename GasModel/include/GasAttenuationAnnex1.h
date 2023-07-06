#ifndef GAS_ATTENUATION_ANNEX1_H
#define GAS_ATTENUATION_ANNEX1_H

#include "Enumerations.h"
#include "GeodeticCoord.h"

class GasAttenuationAnnex1 {
public:
    /// <summary>
    /// Construct an object to calculate gas-related attenuation
    /// </summary>
    /// <param name="ituData">Pointer to object holding all ITU data stores</param>
    GasAttenuationAnnex1() {}
    
    /// <summary>
    /// Section 2.1 "Terrestrial paths"
    /// </summary>
    /// <param name="location">Transmitting location (lat/lon)</param>
    /// <param name="freq_GHz">Transmitting frequency (GHz)</param>
    /// <param name="avail">Desired availability (0-1)</param>
    /// <param name="pathLengthKm">Separation distance between transmitter and receiver (km)</param>
    /// <returns>Attenuation along terrestrial path between transmitter and receiver due to gases (dB)</returns>
    double calculateTerrestrialPathAttenuation_dB(const GeodeticCoord& txLocation, 
                const double& freq_GHz, const double& pathLength_km, 
                const double& rho0_gm3 = 7.5, const bool& useAnnualStandardAtmosphere = false, 
                const Enumerations::Season& season = Enumerations::Season::SummerTime);

    /// <summary>
    /// Section 2.2.1 "Non-negative apparent elevation angles"
    /// </summary>
    /// <param name="txLocation">Transmitting location (lat/lon)</param>
    /// <param name="freq_GHz">Transmitting frequency (GHz)</param>
    /// <param name="avail">Desired availability (0-1)</param>
    /// <param name="elevationAngle_deg">Apparent elevation angle from the transmitting location to the space station (deg)</param>
    /// <returns>Attenuation along the slant path connecting the transmitting Earth station to the receiving space station (dB)</returns>
    double calculateEarthSpaceSlantPathAttenuation_dB(const GeodeticCoord& txLocation, 
                const double& freq_GHz, const double& elevationAngle_deg, 
                const double& rho0_gm3 = 7.5, const bool& useAnnualStandardAtmosphere = false, 
                const Enumerations::Season& season = Enumerations::Season::SummerTime);
};

#endif // GAS_ATTENUATION_ANNEX1_H