#ifndef GAS_ATTENUATION_HELPERS_H
#define GAS_ATTENUATION_HELPERS_H

#include "PowerUnitConversionHelpers.h"
#include "GeodeticCoord.h"
#include "Enumerations.h"
#include "DataStructures.h"

namespace GasAttenuationHelpers {
    /// <summary>
    /// Simple conversion utility for converting ground-level water vapor density to air pressure due to water vapor (used in all of the reference atmosphere models below)
    /// </summary>
    /// <param name="rho_gm3">Ground-level water vapor density (g/m^3)</param>
    /// <param name="temp_K">Temperature at surface (K)</param>
    /// <returns>Air pressure due to water vapor (hPa)</returns>
    double convertWaterVaporGM3toHPA(const double &rho_gm3, const double &temp_K);

    double calculateImaginaryRefractivity_Oxygen(const double& freq_GHz, 
                const double& dryPressure_hPa, const double& waterVapor_hPa, const double& theta);

	double calculateImaginaryRefractivity_Water(const double& freq_GHz, 
                const double& dryPressure_hPa, const double& waterVapor_hPa, const double& theta);

	/// <summary>
	/// Implements the mean annual global reference atmosphere, which approximates the US Standard Atmosphere (1976)
	/// </summary>
	/// <param name="height_km">Height above ground level (km)</param>
	/// <param name="temp_K">Temperature at surface (K)</param>
	/// <param name="dryPressure_hPa">Air pressure not due to water vapor (hPa)</param>
	/// <param name="waterVapor_hPa">Air pressure due to water vapor (hPa)</param>
	/// <param name="rho0_gm3">Ground-level water vapor density (g/m^3)</param>
    void setAtmosphericTermsForUsLocation(const double &height_km, 
                double& temp_K, double& totalPressure_hPa, double& waterVapor_hPa, const double &rho0_gm3 = 7.5);

    /// <summary>
    /// Implements Section 2 of Annex 1 for "Low-latitude annual reference atmosphere"
    /// </summary>
    /// <param name="height_km">Height above ground level (km)</param>
    /// <param name="temp_K">Temperature at surface (K)</param>
    /// <param name="dryPressure_hPa">Air pressure not due to water vapor (hPa)</param>
    /// <param name="rho_gm3">Ground-level water vapor density (g/m^3)</param>
    void setAtmosphericTermsForUsLowLatitude(const double& height_km, 
                double& temp_K, double& totalPressure_hPa, double& rho_gm3);

	/// <summary>
	/// Implements Section 3 of Annex 1 for "Mid-latitude reference atmosphere"
	/// </summary>
     /// <param name="height_km">Height above ground level (km)</param>
    /// <param name="temp_K">Temperature at surface (K)</param>
    /// <param name="dryPressure_hPa">Air pressure not due to water vapor (hPa)</param>
    /// <param name="rho_gm3">Ground-level water vapor density (g/m^3)</param>
	/// <param name="season">Current season (data varies by season)</param>
    void setAtmosphericTermsForUsMidLatitude(const double& height_km, 
                double& temp_K, double& totalPressure_hPa, double& rho_gm3, const Enumerations::Season &season);

    /// <summary>
    /// Implements Section 4 of Annex 1 for "High latitude reference atmosphere"
    /// </summary>
    /// <param name="height_km">Height above ground level (km)</param>
    /// <param name="temp_K">Temperature at surface (K)</param>
    /// <param name="dryPressure_hPa">Air pressure not due to water vapor (hPa)</param>
    /// <param name="rho_gm3">Ground-level water vapor density (g/m^3)</param>
    /// <param name="season">Current season (data varies by season)</param>
    void setAtmosphericTermsForUsHighLatitude(const double& height_km, 
                double& temp_K, double& totalPressure_hPa, double& rho_gm3, const Enumerations::Season& season);
    /// <summary>
    /// Selects the appropriate atmospheric model to use for a given location and season, depending on user preference
    /// </summary>
    /// <param name="location">Location of the transmitting terminal (lat/lon), including its height above ground level (km)</param>
    /// <param name="temp_K">Temperature at surface (K)</param>
    /// <param name="dryPressure_hPa">Air pressure not due to water vapor (hPa)</param>
    /// <param name="waterVapor_hPa">Air pressure due to water vapor (hPa)</param>
    /// <param name="season">Current season (data varies by season)</param>
    /// <param name="useAnnualStandardAtmosphere">Indicates whether the user prefers the mean annual global model or one of the more specific models</param>
    /// <param name="rho0_gm3">Ground-level water vapor density (g/m^3) -- only used when using the mean annual global reference atmosphere</param>
    void setSeasonalAtmosphericTermsForUsLocation(const GeodeticCoord& location, 
                double& temp_K, double& totalPressure_hPa, double& waterVapor_hPa, 
                const Enumerations::Season& season,
                const bool &useAnnualStandardAtmosphere = false, const double &rho0_gm3 = 7.5);

    /// <summary>
    /// Section 1 "specific attenuation" of ITU-R P.676-12
    /// </summary>
    /// <param name="location">Transmitting location (lat/lon)</param>
    /// <param name="freq_GHz">Transmitting frequency (GHz)</param>
    /// <param name="avail">Desired availability (0-1)</param>
    /// <returns>Specific gaseous attenuation due to water at a given transmitting location (dB/km)</returns>
    double calculateSpecificWaterAttenuation_dBPerKm(const double& freq_GHz, 
                const double& temp_K, const double& totalPressure_hPa, const double& waterVapor_hPa);

    /// <summary>
    /// Section 1 "specific attenuation" of ITU-R P.676-12
    /// </summary>
    /// <param name="location">Transmitting location (lat/lon)</param>
    /// <param name="freq_GHz">Transmitting frequency (GHz)</param>
    /// <param name="avail">Desired availability (0-1)</param>
    /// <returns>Specific gaseous attenuation due to oxygen at a given transmitting location (dB/km)</returns>
    double calculateSpecificOxygenAttenuation_dBPerKm(const double& freq_GHz, 
                const double& temp_K, const double& totalPressure_hPa, const double& waterVapor_hPa);

    /// <summary>
    /// Section 1 "specific attenuation" of ITU-R P.676-12
    /// </summary>
    /// <param name="location">Transmitting location (lat/lon)</param>
    /// <param name="freq_GHz">Transmitting frequency (GHz)</param>
    /// <param name="avail">Desired availability (0-1)</param>
    /// <returns>Specific gaseous attenuation at a given transmitting location (dB/km)</returns>
    double calculateSpecificTotalAttenuation_dBPerKm(const double& freq_GHz, 
                const double& temp_K, const double& totalPressure_hPa, const double& waterVapor_hPa);
	/// <summary>
	/// Implements radio refractive index formula as described in ITU-R P.453-14
	/// </summary>
	/// <param name="dryPressure_hPa"></param>
	/// <param name="waterVaporPressure_hPa"></param>
	/// <param name="temp_K"></param>
	/// <returns></returns>
    double calculateRadioRefractiveIndex(const double &totalPressure_hPa, 
                const double &waterVaporPressure_hPa, const double &temp_K);

    /* ITU-R P.676-12 ANNEX 2 HELPERS START HERE */
    double calculateHeightOfWaterComponent_km(const double& freq_GHz, 
                const double& temp_K, const double& rho0_gm3, const double& pressureRatio);

    double calculateHeightOfOxygenComponent_km(const double& freq_GHz, 
                const double& temp_K, const double& pressureRatio);

    /// <summary>
    /// Implementation of Section 2.3 "Zenith path water vapour attenuation" from ITU-R P.676-12
    /// </summary>
    /// <param name="freq_GHz"></param>
    /// <param name="integratedWaterVapor_kgm2"></param>
    /// <returns></returns>
    double calculateZenithWaterVaporAttenuation_dB(const double &freq_GHz, 
                const double &height_km, const double &integratedWaterVapor_kgm2);

} // end namespace GasAttenuationHelpers

#endif /* GAS_ATTENUATION_HELPERS_H */