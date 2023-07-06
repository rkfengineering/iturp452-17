#include "GasAttenuationAnnex2.h"
#include "GasAttenuationHelpers.h"

#include <iostream>
#include <limits>

double GasAttenuationAnnex2::calculateZenithAttenuation_dB(const GeodeticCoord& txLocation, const double& freq_GHz, 
                const double& rho0_gm3, const Enumerations::Season& season) {
	if (txLocation.height_km < 0.0 || txLocation.height_km > 10.0) {
		std::ostringstream oStrStream;
		oStrStream << "GasAttenuationAnnex2::calculateZenithAttenuation_dB(): " 
					<< "The given height is outside of the valid range for heights [0, 10] km: " 
					<< txLocation.height_km << " km!";
		throw std::domain_error(oStrStream.str());
	}

	// Find atmospheric parameters for the current location
	double temp_K, totalPressure_hPa, waterVapor_hPa;
	GasAttenuationHelpers::setSeasonalAtmosphericTermsForUsLocation(txLocation, temp_K, totalPressure_hPa, waterVapor_hPa, season, false, rho0_gm3);

	// Total pressure is the super of the dry pressure and water vapor pressure
	const double PRESSURE_RATIO = totalPressure_hPa / 1013.25;

	const double height_oxygenComp_km = GasAttenuationHelpers::calculateHeightOfOxygenComponent_km(freq_GHz, temp_K, PRESSURE_RATIO);

	const double height_waterComp_km = GasAttenuationHelpers::calculateHeightOfWaterComponent_km(freq_GHz, temp_K, rho0_gm3, PRESSURE_RATIO);

	const double waterSpecificAttenuation_dBPerKm = GasAttenuationHelpers::calculateSpecificWaterAttenuation_dBPerKm(freq_GHz, temp_K, totalPressure_hPa, waterVapor_hPa);
	const double oxygenSpecificAttenuation_dBPerKm = GasAttenuationHelpers::calculateSpecificOxygenAttenuation_dBPerKm(freq_GHz, temp_K, totalPressure_hPa, waterVapor_hPa);

	// Equation #39
	const double ZENITH_ATTENUATION_DB = oxygenSpecificAttenuation_dBPerKm * height_oxygenComp_km + waterSpecificAttenuation_dBPerKm * height_waterComp_km;

	return ZENITH_ATTENUATION_DB;
}

double GasAttenuationAnnex2::calculateEarthSpacePathAttenuation_dB(const double& height_km, const double& freq_GHz, 
                const double& elevationAngle_deg, const double& temp_K, const double& totalPressure_hPa, 
                const double& waterPressure_hPa, const double& integratedWaterVapor_kgm2) {
	if (elevationAngle_deg < 5.0 || elevationAngle_deg > 90.0) {
		std::ostringstream oStrStream;
		oStrStream << "GasAttenuationAnnex2::calculateEarthSpacePathAttenuation_dB(): The given elevation angle is outside of the valid range [5, 90] deg: " << elevationAngle_deg << " deg!";
		throw std::domain_error(oStrStream.str());
	}

	// Total pressure is the super of the dry pressure and water vapor pressure
	const double PRESSURE_RATIO = totalPressure_hPa / 1013.25;

	const double height_oxygenComp_km = GasAttenuationHelpers::calculateHeightOfOxygenComponent_km(freq_GHz, temp_K, PRESSURE_RATIO);
	const double oxygenSpecificAttenuation_dBPerKm = GasAttenuationHelpers::calculateSpecificOxygenAttenuation_dBPerKm(freq_GHz, temp_K, totalPressure_hPa, waterPressure_hPa);

	const double waterVaporAttenuation_dB = GasAttenuationHelpers::calculateZenithWaterVaporAttenuation_dB(freq_GHz, height_km, integratedWaterVapor_kgm2);

	// Equation #41
	const double TOTAL_ATTEN_DB = (height_oxygenComp_km * oxygenSpecificAttenuation_dBPerKm + waterVaporAttenuation_dB) / std::sin(MathHelpers::deg2rad(elevationAngle_deg));

	return TOTAL_ATTEN_DB;
}

double GasAttenuationAnnex2::calculateEarthSpacePathAttenuationSeasonal_dB(const GeodeticCoord& txLocation, 
                const double& freq_GHz, const double& rho0_gm3, const double& integratedWaterVapor_kgm2, 
                const double& elevationAngle_deg, const Enumerations::Season& season) {
	// Calculate atmospheric parameters for the current season
	double temp_K, totalPressure_hPa, waterVapor_hPa;
	GasAttenuationHelpers::setSeasonalAtmosphericTermsForUsLocation(txLocation, temp_K, totalPressure_hPa, waterVapor_hPa, season, false, rho0_gm3);

	const double TOTAL_ATTEN_DB = calculateEarthSpacePathAttenuation_dB(txLocation.height_km, freq_GHz, 
				elevationAngle_deg, temp_K, 
				totalPressure_hPa, waterVapor_hPa, integratedWaterVapor_kgm2);
	return TOTAL_ATTEN_DB;
}

double GasAttenuationAnnex2::calculateEarthInclinedPathAttenuation_dB(const GeodeticCoord& txLocation, 
                const GeodeticCoord& rxLocation, const double& freq_GHz, const double& rho0_gm3, 
                const double& elevationAngle_deg, const Enumerations::Season& season) {
	if (elevationAngle_deg < 5.0 || elevationAngle_deg > 90.0) {
		std::ostringstream oStrStream;
		oStrStream << "GasAttenuationAnnex2::calculateEarthInclinedPathAttenuation_dB(): " 
					<< "The given elevation angle is outside of the valid range [5, 90] deg: " 
					<< elevationAngle_deg << " deg!";
		throw std::domain_error(oStrStream.str());
	}
	if (txLocation.height_km < 0.0 || txLocation.height_km > 10.0) {
		std::ostringstream oStrStream;
		oStrStream << "GasAttenuationAnnex2::calculateEarthInclinedPathAttenuation_dB(): " 
					<< "The given transmitter height is outside of the valid range for heights [0, 10] km: " 
					<< txLocation.height_km << " km!";
		throw std::domain_error(oStrStream.str());
	}
	if (rxLocation.height_km < 0.0 || rxLocation.height_km > 10.0) {
		std::ostringstream oStrStream;
		oStrStream << "GasAttenuationAnnex2::calculateEarthInclinedPathAttenuation_dB(): " 
					<< "The given receiver height is outside of the valid range for heights [0, 10] km: " 
					<< rxLocation.height_km << " km!";
		throw std::domain_error(oStrStream.str());
	}

	// Attempt to find a rho0_gm3 for the current location (and current height)
	double temp_K, totalPressure_hPa, waterVapor_hPa;
	GasAttenuationHelpers::setSeasonalAtmosphericTermsForUsLocation(txLocation, temp_K, totalPressure_hPa, waterVapor_hPa, season, false, rho0_gm3);

	// Total pressure is the super of the dry pressure and water vapor pressure
	const double PRESSURE_RATIO = totalPressure_hPa / 1013.25;

	const double height_oxygenComp_km = GasAttenuationHelpers::calculateHeightOfOxygenComponent_km(freq_GHz, temp_K, PRESSURE_RATIO);
	// Equation #42
	const double TX_HEIGHT_KM = txLocation.height_km;
	const double RX_HEIGHT_KM = rxLocation.height_km;
	const double height_prime_oxygenComp_km = height_oxygenComp_km * (exp(-TX_HEIGHT_KM / height_oxygenComp_km) - exp(-RX_HEIGHT_KM / height_oxygenComp_km));

	const double height_waterComp_km = GasAttenuationHelpers::calculateHeightOfWaterComponent_km(freq_GHz, temp_K, rho0_gm3, PRESSURE_RATIO);
	// Equation #43
	const double height_prime_waterComp_km = height_waterComp_km * (exp(-TX_HEIGHT_KM / height_waterComp_km) - exp(-RX_HEIGHT_KM / height_waterComp_km));

	const double waterSpecificAttenuation_dBPerKm = GasAttenuationHelpers::calculateSpecificWaterAttenuation_dBPerKm(freq_GHz, temp_K, totalPressure_hPa, waterVapor_hPa);
	const double oxygenSpecificAttenuation_dBPerKm = GasAttenuationHelpers::calculateSpecificOxygenAttenuation_dBPerKm(freq_GHz, temp_K, totalPressure_hPa, waterVapor_hPa);

	// Equation #39
	const double INCLINED_ATTENUATION_DB = oxygenSpecificAttenuation_dBPerKm * height_prime_oxygenComp_km + waterSpecificAttenuation_dBPerKm * height_prime_waterComp_km;

	return INCLINED_ATTENUATION_DB;
}