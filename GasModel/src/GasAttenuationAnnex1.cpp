#include "GasAttenuationAnnex1.h"
#include "GasAttenuationHelpers.h"
#include "MathHelpers.h"
#include "PhysicalConstants.h"

#include <cmath>
#include <iostream>
#include <ostream>
#include <sstream>

double GasAttenuationAnnex1::calculateTerrestrialPathAttenuation_dB(const GeodeticCoord& txLocation, 
                const double& freq_GHz, const double& pathLength_km, 
                const double& rho0_gm3, const bool& useAnnualStandardAtmosphere, 
                const Enumerations::Season& season) {
	
	// Require seasonal input when calculating atmospheric terms
	double temp_K, totalPressure_hPa, waterVapor_hPa;
	GasAttenuationHelpers::setSeasonalAtmosphericTermsForUsLocation(txLocation, temp_K, totalPressure_hPa, waterVapor_hPa, season, useAnnualStandardAtmosphere, rho0_gm3);
	
	const double specificAttenuation_dBPerKm = GasAttenuationHelpers::calculateSpecificTotalAttenuation_dBPerKm(freq_GHz, temp_K, totalPressure_hPa, waterVapor_hPa);
	// Equation #10 from ITU-R P.676-12
	const double terrestialPathLoss_dB = specificAttenuation_dBPerKm * pathLength_km;

	return terrestialPathLoss_dB;
}

double GasAttenuationAnnex1::calculateEarthSpaceSlantPathAttenuation_dB(const GeodeticCoord& txLocation, 
                const double& freq_GHz, const double& elevationAngle_deg, 
                const double& rho0_gm3, const bool& useAnnualStandardAtmosphere, 
                const Enumerations::Season& season) {
	if (txLocation.height_km < 0.0 || txLocation.height_km > 100.0) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: GasAttenuationAnnex1::calculateEarthSpaceSlantPathAttenuation_dB(): " << 
			"The given height falls outside of the range of heights ([0-100] km) specified in the ITU's standard atmospheric profiles: " 
			<< txLocation.height_km << " km!" << std::endl;
		throw std::domain_error(oStrStream.str());
	}
	
	// If one of the six reference standard atmospheres specified in ITU-R P.835 is used, then the atmospheric profile is defined for geometric heights up to 100 km
	// This means that the maximum number of layers is 922
	// NOTE: If we use other (local) data, then we need to calculate these terms ourselves
	const uint64_t MAX_LAYER_IND = 922;

	const double EXP_SINGLE_LAYER_FACTOR = std::exp(0.01);
	// Equation #16a from ITU-R P.676-12
	const uint64_t MIN_LAYER_IND = static_cast<uint64_t>(std::floor(100.0 * std::log(10.0e4 * txLocation.height_km * (EXP_SINGLE_LAYER_FACTOR - 1.0) + 1.0) + 1.0));

	if (MAX_LAYER_IND - MIN_LAYER_IND < 50) {
		/* TODO: Setup a verbose option to allow warnings to be displayed
		std::cerr << "WARNING: GasAttenuationAnnex1::calculateEarthSpaceSlantPathAttenuation_dB(): " 
					<< "There is possible degraded accuracy for slant paths with small altitude changes " 
					<< "(i.e., between two airborne platforms)" << std::endl;
		*/
	}

	// Create a copy of the Tx location so that we can increment its height as we step through the layers
	double layerHeight_km = 0.0;
	// Store atmospheric conditions and other data for each layer
	std::vector<double> layerThicknessKmList;
	layerThicknessKmList.reserve(MAX_LAYER_IND);
	std::vector<double> layerRefractIndList;
	layerRefractIndList.reserve(MAX_LAYER_IND);
	std::vector<double> specificAttenuationList;
	specificAttenuationList.reserve(MAX_LAYER_IND);
	double layerThickness_km, layerMidpointHeight_km, layerRefractInd, specificAttenuation_dBPerKm;
	double temp_K, totalPressure_hPa, waterVapor_hPa;
	// Starting layer index at 1, so shifting MAX_LAYER_IND up by 1
	for (uint64_t layerInd = 1; layerInd < MAX_LAYER_IND + 1; layerInd++) {
		// Equation #14 from ITU-R P.676-12
		layerThickness_km = 0.0001 * std::exp((layerInd - 1) / 100.0);
		layerThicknessKmList.push_back(layerThickness_km);
		
		layerMidpointHeight_km = layerHeight_km + layerThickness_km / 2.0;
		GeodeticCoord layerCoord(txLocation.m_longitude_deg, txLocation.m_latitude_deg, layerMidpointHeight_km);
		// Attempt to find a rho0_gm3 for the current location (and current height)
		GasAttenuationHelpers::setSeasonalAtmosphericTermsForUsLocation(layerCoord, temp_K, totalPressure_hPa, waterVapor_hPa, 
			season, useAnnualStandardAtmosphere, rho0_gm3);

		// Using these atmospheric terms, calculate the refractive index associated with the lowest layer
		layerRefractInd = GasAttenuationHelpers::calculateRadioRefractiveIndex(totalPressure_hPa, waterVapor_hPa, temp_K);
		layerRefractIndList.push_back(layerRefractInd);

		specificAttenuation_dBPerKm = GasAttenuationHelpers::calculateSpecificTotalAttenuation_dBPerKm(freq_GHz,
			temp_K, totalPressure_hPa, waterVapor_hPa);
		specificAttenuationList.push_back(specificAttenuation_dBPerKm);

		// Increment the layer's height to the next layer
		layerHeight_km += layerThickness_km;
	}

	// Defined as Beta_1 used in Equation #19b from ITU-R P.676-12
	const double COMPLEMENT_ELEVATION_RAD = MathHelpers::deg2rad(90.0 - elevationAngle_deg);
	
	double layerZenithAngle = COMPLEMENT_ELEVATION_RAD;
	double pathLength_km, layerAlphaAngle, tempInput;
	double totalAttenuation_dB = 0.0;
	double currentLayerHeightFromEarthCenter_km = PhysicalConstants::EARTH_MEAN_RADIUS_KM;
	for (uint64_t layerInd = 0; layerInd < MAX_LAYER_IND - 1; layerInd++) {
		while ((layerInd + 1) < MIN_LAYER_IND) {
			currentLayerHeightFromEarthCenter_km += layerThicknessKmList[layerInd];
			layerInd++;
		}

		// Equation #18b from ITU-R P.676-12
		tempInput = (currentLayerHeightFromEarthCenter_km * std::sin(layerZenithAngle)) / (currentLayerHeightFromEarthCenter_km + layerThicknessKmList[layerInd]);
		// std::asin is only defined for [-1,1]
		layerAlphaAngle = std::asin(std::min(1.0, tempInput));

		// Equation #17 from ITU-R P.676-12
		tempInput = std::sqrt(MathHelpers::simpleSquare(currentLayerHeightFromEarthCenter_km * std::cos(layerZenithAngle)) +
			2.0 * currentLayerHeightFromEarthCenter_km * layerThicknessKmList[layerInd] + MathHelpers::simpleSquare(layerThicknessKmList[layerInd]));
		pathLength_km = -currentLayerHeightFromEarthCenter_km * std::cos(layerZenithAngle) + tempInput;

		// Update this angle to the next layer's zenith angle -- Equation #19a from ITU-R P.676-12
		tempInput = (layerRefractIndList[layerInd] / layerRefractIndList[layerInd + 1]) * std::sin(layerAlphaAngle);
		// std::asin is only defined for [-1,1]
		layerZenithAngle = std::asin(std::min(1.0, tempInput));

		// Equation #13 from ITU-R P.676-12
		totalAttenuation_dB += specificAttenuationList[layerInd] * pathLength_km;

		// Increment to the next layer and update to that layer's height
		currentLayerHeightFromEarthCenter_km += layerThicknessKmList[layerInd];
	}
	
	return totalAttenuation_dB;
}