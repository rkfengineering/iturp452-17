#include "DataStructures.h"
#include "GasAttenuation.h"
#include "GasAttenuationAnnex1.h"
#include "GasAttenuationAnnex2.h"
#include "GasAttenuationHelpers.h"

double GasAttenuation::performCalculations(const Gas::DataLoader& dataLoader, 
            const GeodeticCoord& location, const double& freq_GHz, const double& elevationAngle_deg, 
            const double& exceedance, const Enumerations::Season& season,
            const bool& useAnnexOne, const bool& useStandardAtmosphere) {
    // Equations #61-62 from ITU-R P.618-13
	const double RHO_GM3 = dataLoader.fetchSurfaceWaterVaporDensity(location, exceedance);

	if (useAnnexOne) {
		GasAttenuationAnnex1 gasAttenuationAnnex1;
		// Perform calculation
		return gasAttenuationAnnex1.calculateEarthSpaceSlantPathAttenuation_dB(location, freq_GHz, 
			elevationAngle_deg, RHO_GM3, useStandardAtmosphere, season);
	}
	else {
		GasAttenuationAnnex2 gasAttenuationAnnex2;
		// Set up input data for atmospheric conditions at this location/exceedance
		double temp_K, totalPressure_hPa, wetPressure_hPa;
		GasAttenuationHelpers::setAtmosphericTermsForUsLocation(location.height_km, temp_K, totalPressure_hPa, wetPressure_hPa, RHO_GM3);
		const double INTEGRATED_WATER_VAPOR_CONTENT = dataLoader.fetchTotalWaterVaporContent(location, exceedance);
		
		// Perform calculation
		return gasAttenuationAnnex2.calculateEarthSpacePathAttenuation_dB(location.height_km, freq_GHz,
			 elevationAngle_deg, temp_K, totalPressure_hPa, wetPressure_hPa, INTEGRATED_WATER_VAPOR_CONTENT);
	}
}