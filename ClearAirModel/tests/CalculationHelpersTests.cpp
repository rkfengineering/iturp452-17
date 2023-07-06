#include "gtest/gtest.h"

#include "GasAttenuationHelpers.h"

namespace {
	// Use when expected an exact match
	double constexpr TOLERANCE = 1.0e-6;
}

TEST(GasAttenuationHelpersTests, GasAttenuation_imaginaryRefractivityOxygenTest) {
	const std::vector<double> FREQ_LIST = {
		12.0, 20.0, 60.0, 90.0, 130.0
	};
	const double DRY_PRESSURE_HPA = 1013.25;
	const double WATER_PRESSURE_HPA = 9.97288879;
	const double TEMP_K = 288.15;
	const double THETA = 300.0 / TEMP_K;

	const std::vector<double> EXPECTED_N_DOUBLE_PRIME_LIST = {
		0.0039827, 0.00326471, 1.33914604, 0.00237300, 0.00175440
	};

	for (uint16_t freqInd = 0; freqInd < FREQ_LIST.size(); freqInd++) {
		const double N_DOUBLE_PRIME = GasAttenuationHelpers::calculateImaginaryRefractivity_Oxygen(FREQ_LIST[freqInd],
			DRY_PRESSURE_HPA, WATER_PRESSURE_HPA, THETA);
		EXPECT_NEAR(EXPECTED_N_DOUBLE_PRIME_LIST[freqInd], N_DOUBLE_PRIME, TOLERANCE);
	}
}

TEST(GasAttenuationHelpersTests, GasAttenuation_imaginaryRefractivityWaterTest) {
	const std::vector<double> FREQ_LIST = {
		12.0, 20.0, 60.0, 90.0, 130.0
	};
	const double DRY_PRESSURE_HPA = 1013.25;
	const double WATER_PRESSURE_HPA = 9.97288879;
	const double TEMP_K = 288.15;
	const double THETA = 300.0 / TEMP_K;

	const std::vector<double> EXPECTED_N_DOUBLE_PRIME_LIST = {
		0.00436602, 0.02666135, 0.01417966, 0.02087750, 0.03177704
	};

	for (uint16_t freqInd = 0; freqInd < FREQ_LIST.size(); freqInd++) {
		const double N_DOUBLE_PRIME = GasAttenuationHelpers::calculateImaginaryRefractivity_Water(FREQ_LIST[freqInd],
			DRY_PRESSURE_HPA, WATER_PRESSURE_HPA, THETA);
		EXPECT_NEAR(EXPECTED_N_DOUBLE_PRIME_LIST[freqInd], N_DOUBLE_PRIME, TOLERANCE);
	}
}