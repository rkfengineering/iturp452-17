#include "gtest/gtest.h"

#include "AttenuationTests.h"
#include "GasAttenuation.h"
#include "GasAttenuationAnnex1.h"
#include "GasAttenuationAnnex2.h"
#include "GeodeticCoord.h"

namespace {
	// Use when expected an exact match
	double constexpr TOLERANCE = 1.0e-6;
}

TEST_F(AttenuationTests, GasAttenuation_calculateAttenuationTest) {
	// From ITU validation spreadsheet titled "CG-3M3J-13-ValEx-Rev6.1.2, page "P618-13 A_Gas_A1_2.2.1a"
	// Arrange
	const double FREQ_GHZ = 28.0;
	const double HEIGHT_KM = 0.0;
	const double ELEV_ANGLE_DEG = 30.0;
	const double RHO_GM3 = 7.5;
	const GeodeticCoord COORD(-80.0, 38.0, HEIGHT_KM);
	
	const double EXPECTED_GAS_ATTEN = 0.47081174;

	// Act
	GasAttenuationAnnex1 gasAttenuation;
	const double GAS_ATTEN_DB = gasAttenuation.calculateEarthSpaceSlantPathAttenuation_dB(COORD, FREQ_GHZ,
		ELEV_ANGLE_DEG, RHO_GM3, true);

	// Assert
	EXPECT_NEAR(EXPECTED_GAS_ATTEN, GAS_ATTEN_DB, TOLERANCE);
}

TEST_F(AttenuationTests, GasAttenuation_calculateSummerAttenuationTest) {
	// Arrange
	const double longitude_deg = -85.0;
	const double latitude_deg = 30.0;
	const double height_km = 0.0;
	const double freq_GHz = 17.3;
	const double elevationAngle_deg = 20.0;
	const GeodeticCoord location(longitude_deg, latitude_deg, height_km);

	// NOTE: This expected value was copied from the output observed, rather than from an ITU validation source
	const double EXPECTED_GAS_ATTEN_DB = 0.5196453638;

	// Act
	GasAttenuationAnnex1 gasAttenuation;
	const double ATTEN_DB = gasAttenuation.calculateEarthSpaceSlantPathAttenuation_dB(location, freq_GHz, elevationAngle_deg, 
		7.5, false, Enumerations::Season::SummerTime);

	// Assert
	EXPECT_NEAR(EXPECTED_GAS_ATTEN_DB, ATTEN_DB, TOLERANCE);
}

TEST_F(AttenuationTests, GasAttenuation_calculateWinterAttenuationTest) {
	// Arrange
	const double longitude_deg = -85.0;
	const double latitude_deg = 30.0;
	const double height_km = 0.0;
	const double freq_GHz = 17.3;
	const double elevationAngle_deg = 20.0;
	const GeodeticCoord location(longitude_deg, latitude_deg, height_km);

	// NOTE: This expected value was copied from the output observed, rather than from an ITU validation source
	const double EXPECTED_GAS_ATTEN_DB = 0.2714642853;

	// Act
	GasAttenuationAnnex1 gasAttenuation;
	const double ATTEN_DB = gasAttenuation.calculateEarthSpaceSlantPathAttenuation_dB(location, freq_GHz, elevationAngle_deg, 
		7.5, false, Enumerations::Season::WinterTime);

	// Assert
	EXPECT_NEAR(EXPECTED_GAS_ATTEN_DB, ATTEN_DB, TOLERANCE);
}

TEST_F(AttenuationTests, GasAttenuationAnnex2_basicTest) {
	// From ITU validation spreadsheet titled "CG-3M3J-13-ValEx-Rev6.1.2, page "P676-12 A_Gas_A2"
	// Arrange
	const std::vector<double> HEIGHT_LIST = { 
		0.03138298, 0.04612299, 0.00000000, 0.03138298, 0.04612299, 0.00000000, 0.03138298, 0.04612299, 
		0.00000000, 0.03138298, 0.04612299, 0.00000000, 0.03138298, 0.04612299, 0.00000000, 0.03138298, 
		0.04612299, 0.00000000, 0.03138298, 0.04612299, 0.00000000
	};
	const std::vector<double> ELEVATION_LIST = {
		31.076991236,40.232035996,46.359692612,31.076991236,40.232035996,46.359692612,31.076991236,
		40.232035996,46.359692612,31.076991236,40.232035996,46.359692612,31.076991236,40.232035996,
		46.359692612,31.076991236,40.232035996,46.359692612,31.076991236,40.232035996,46.359692612
	};
	const std::vector<double> FREQ_LIST = { 14.250000000,14.250000000,14.250000000,
		14.250000000,14.250000000,14.250000000,14.250000000,14.250000000,
		14.250000000,14.250000000,14.250000000,14.250000000,29.000000000,29.000000000,29.000000000,
		29.000000000,29.000000000,29.000000000,29.000000000,29.000000000,29.000000000
	};
	const std::vector<double> TEMP_LIST = { 283.610875556,
		288.089736889,293.369679467,283.610875556,288.089736889,293.369679467,283.610875556,288.089736889,293.369679467,
		283.610875556,288.089736889,293.369679467,283.610875556,288.089736889,293.369679467,283.610875556,
		288.089736889,293.369679467,283.610875556,288.089736889,293.369679467
	};
	const std::vector<double> WATER_PRESSURE_LIST = {
		18.056519975,24.278798976,30.772004321,18.862715395,24.925956412,31.407615120,19.348351519,
		25.382396240,31.834258792,19.685301134,25.720443531,32.165535028,18.056519975,24.278798976,
		30.772004321,18.862715395,24.925956412,31.407615120,19.348351519,25.382396240,31.834258792
	};
	const std::vector<double> DRY_PRESSURE_LIST = {
		1009.485611798,1007.721474334,1013.250000000,1009.485611798,1007.721474334,1013.250000000,1009.485611798,1007.721474334,1013.250000000,
		1009.485611798,1007.721474334,1013.250000000,1009.485611798,1007.721474334,
		1013.250000000,1009.485611798,1007.721474334,1013.250000000,1009.485611798,1007.721474334,1013.250000000
	};
	const std::vector<double> VAPOR_CONTENT_LIST = {
		33.729465269,36.04810935,37.95559991,35.72555437,37.51301422,40.51190587,37.29595324,38.32964667,
		42.21022387,38.24407811,39.05128132,43.41046058,33.72946527,36.04810935,
		37.95559991,35.72555437,37.51301422,40.51190587,37.29595324,38.32964667,42.21022387
	};

	const std::vector<double> EXPECTED_ATTEN_LIST = {
		0.22687404,0.18948086,0.17601146,0.23608491,0.19496915,0.18467591,0.24341143,0.19806196,0.19052114,0.24788012,
		0.20080484,0.19469701,0.83765994,0.70696616,0.66597767,0.88047255,0.73238144,0.70648507,0.91453438, 0.74667715, 0.73378408
	};

	// Act
	GasAttenuationAnnex2 gasAttenuation;
	for (uint16_t runInd = 0; runInd < HEIGHT_LIST.size(); runInd++) {
		const double totalPressure_hPa = DRY_PRESSURE_LIST[runInd] + WATER_PRESSURE_LIST[runInd];
		const double ACTUAL_ATTEN_DB = gasAttenuation.calculateEarthSpacePathAttenuation_dB(HEIGHT_LIST[runInd], FREQ_LIST[runInd],
			ELEVATION_LIST[runInd], TEMP_LIST[runInd], totalPressure_hPa, WATER_PRESSURE_LIST[runInd], VAPOR_CONTENT_LIST[runInd]);

		// Assert
		EXPECT_NEAR(EXPECTED_ATTEN_LIST[runInd], ACTUAL_ATTEN_DB, TOLERANCE);
	}
}

TEST_F(AttenuationTests, FetchSurfaceWaterVaporDensityTest) {
	// From ITU validation spreadsheet titled "CG-3M3J-13-ValEx-Rev6.1.2, page "P836-6 WV"
	// Arrange
	const std::vector<double> LAT_LIST = {
		3.133, 3.133, 3.133, 3.133, 22.9, 22.9, 22.9, 22.9, 23.00, 
		23.00, 23.00, 23.00, 25.78, 25.78, 25.78, 25.78, 28.717, 
		28.717, 28.717, 28.717
	};
	const std::vector<double> LON_LIST = {
		101.7, 101.7, 101.7, 101.7, -43.23, -43.23, -43.23, -43.23, 
		30.00, 30.00, 30.00, 30.00, -80.22, -80.22, -80.22, -80.22, 
		77.3, 77.3, 77.3, 77.3
	};
	const std::vector<double> HEIGHT_METERS_LIST = {
		51.25145595, 51.25145595, 51.25145595, 51.25145595, 0.00000000,
		0.00000000, 0.00000000, 0.00000000, 187.59375000, 187.59375000,
		187.59375000, 187.59375000, 8.61728000, 8.61728000, 8.61728000,
		8.61728000, 209.38369895, 209.38369895, 209.38369895, 209.38369895
	};
	const std::vector<double> PERCENT_EXCEED_LIST = {
		0.1, 0.15, 0.3, 0.35, 0.1, 0.15, 0.3, 0.35,
		0.1, 0.15, 0.3, 0.35, 0.1, 0.15, 0.3, 0.35,
		0.1, 0.15, 0.3, 0.35
	};
	const std::vector<double> EXPECTED_WATER_VAPOR_DENSITY_LIST = {
		24.32302408, 24.19572903, 23.95244058, 23.89334611, 21.59164912, 
		21.46164369, 21.24753319, 21.18676013, 11.98314750, 11.72133345, 
		11.22988725, 11.11750130, 23.44828894, 23.28196028, 23.00080897, 
		22.93780875, 26.00984124, 25.76970603, 25.39894305, 25.31476962
	};

	// Act
	const Gas::DataLoader dataLoader(*m_ituDataFileReader);
	for (uint32_t latInd = 0; latInd < LAT_LIST.size(); latInd++) {
		const GeodeticCoord CURRENT_COORD(LON_LIST[latInd], LAT_LIST[latInd], HEIGHT_METERS_LIST[latInd] / 1.0e3);
		const double WATER_VAPOR_DENSITY = dataLoader.fetchSurfaceWaterVaporDensity(CURRENT_COORD, 
			PERCENT_EXCEED_LIST[latInd] / 100.0);

		// Assert
		EXPECT_NEAR(EXPECTED_WATER_VAPOR_DENSITY_LIST[latInd], WATER_VAPOR_DENSITY, TOLERANCE);
	}
}