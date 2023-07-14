#include "gtest/gtest.h"
#include "PathProfile.h"
#include "InvCumNorm.h"
#include "EffectiveEarth.h"
#include <filesystem>

//Example Profile Path from ITU validation spreadsheet titled "delB_valid_temp.xlsx", pages "Path 1" to "Path 4"
//embedded in ITU validation document titled "Validation Examples for the delta Bullington diffraction prediction method"

//Unit Tests for smaller functions (not using ITU Validation data)

const std::filesystem::path testPathFileDir("/home/ayeh/itu/ituModels/iturp452/ClearAirModel/tests/test_paths");

namespace {
	// Use when expected an exact match
	double constexpr TOLERANCE = 1.0e-6;
}

TEST(ProfilePathTests, ProfilePath_loadDHProfileTest){
	const std::vector<std::string> PROFILE_LIST = {
		"path1.csv"
	};
	const std::vector<double> EXPECTED_FIRST_D_LIST = {
		0
	};
	const std::vector<double> EXPECTED_FIRST_H_LIST = {
		800
	};
	const std::vector<double> EXPECTED_LAST_D_LIST = {
		83.07993853
	};
	const std::vector<double> EXPECTED_LAST_H_LIST = {
		302
	};
	const std::vector<double> EXPECTED_LENGTH_LIST = {
		1663
	};
	for (uint16_t profileInd = 0; profileInd < PROFILE_LIST.size(); profileInd++) {
		const PathProfile::Path PROFILE(testPathFileDir/std::filesystem::path(PROFILE_LIST[profileInd]));
		EXPECT_NEAR(EXPECTED_FIRST_D_LIST[profileInd], PROFILE.front().d_km, TOLERANCE);
		EXPECT_NEAR(EXPECTED_FIRST_H_LIST[profileInd], PROFILE.front().h_masl, TOLERANCE);
		EXPECT_NEAR(EXPECTED_LAST_D_LIST[profileInd], PROFILE.back().d_km, TOLERANCE);
		EXPECT_NEAR(EXPECTED_LAST_H_LIST[profileInd], PROFILE.back().h_masl, TOLERANCE);
		EXPECT_NEAR(EXPECTED_LENGTH_LIST[profileInd], PROFILE.size(), TOLERANCE);
	}	
}

//Compared against NORM.INV Built-in Macro from Microsoft Excel
TEST(HelpersTests, InvCumNormTest){
	const std::vector<double> INPUTS_LIST = {
		1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5,
		0.01345, 0.42039, 0.0042598, 0.000050938
	};
	const std::vector<double> EXPECTED_RESULTS_LIST = {
		-4.753424309, -4.264890794, -3.719016485,
		-3.090232306, -2.326347874, -1.281551566, 0.0, -2.212965852,
		-0.200895866, -2.630752722, -3.886079871,
	};

	for (uint16_t ind = 0; ind < INPUTS_LIST.size(); ind++) {
		const double OUTPUT = inv_cum_norm(INPUTS_LIST[ind]);
		//Maximum error of 0.00054 from Attachment 3 to Annex 1 ITU-R P.452-17
		EXPECT_NEAR(EXPECTED_RESULTS_LIST[ind], OUTPUT, 0.00054);
	}
	
	//use 1e-6 for all probability below 1e-6
	EXPECT_NEAR(inv_cum_norm(1e-7), inv_cum_norm(1e-6), TOLERANCE);
}

//quick logic test. make sure values make sense (the path isn't extremely jagged)
//Endpoints estimated from Excel Trendline
TEST(EffectiveEarthTests, EffectiveEarthTests_smoothEarthAMSLHeights){
	const double EXPECTED_START = 635.2;
	const double EXPECTED_END = 357.3;

	const PathProfile::Path p(testPathFileDir/std::filesystem::path("path1.csv"));
	const EffectiveEarth::HeightPair EFF_HEIGHT = EffectiveEarth::smoothEarthHeights_AMSL(p);

	//high tolerance since the values were read off a graph
	EXPECT_NEAR(EXPECTED_START,EFF_HEIGHT.tx_val,0.5);
	EXPECT_NEAR(EXPECTED_END,EFF_HEIGHT.rx_val,0.5);
}
//test effective earth radius. check close to 8500 for some DN from data map
