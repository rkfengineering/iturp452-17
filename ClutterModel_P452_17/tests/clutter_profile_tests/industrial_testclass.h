
#ifndef TEST_PROFILE_INDUSTRIAL_H
#define TEST_PROFILE_INDUSTRIAL_H

#include "gtest/gtest.h"

#include <ITUR_P452/PathProfile.h>
#include "ClutterModel_P452_17/ClutterLoss.h"
#include "Common/Enumerations.h"
#include <filesystem>


class IndustrialProfileTests : public testing::Test {
protected:
	IndustrialProfileTests() {}

	// Per-test-suite set-up.
	// Called before the first test in this test suite.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
        K_PATH = PathProfile::Path(CMAKE_CLEARAIR_SRC_DIR 
						/ std::filesystem::path("tests/test_paths")
						/ std::filesystem::path("test_profile_flat_land_5km_Industrial.csv"));
	}

	// Per-test-suite tear-down.
	// Called after the last test in this test suite.
	// Can be omitted if not needed.
	static void TearDownTestCase() { }

	// You can define per-test set-up logic as usual.
	void SetUp() override { }

	// You can define per-test tear-down logic as usual.
	void TearDown() override { }

	// Data management object that should only be initialized once for all unit tests
	//TODO move this stuff to a test case set up class
    static PathProfile::Path K_PATH;

	const std::vector<double> FREQ_GHZ_LIST = {
		2,0.1,0.15,0.225,0.3375,0.50625,0.759375,1.1390625,
		1.70859375,2.562890625,3.844335938,5.766503906,8.649755859,
		12.97463379,19.46195068,29.19292603,43.78938904,50,2,2,2,2,
		2,2,2,2,2,2,2,2,2,2,2,2,2
	};
    const std::vector<double> P_LIST = {
        49,49,49,49,49,49,49,49,49,49,49,49,49,
		49,49,49,49,49,1,4,7,10,13,16,19,22,
		25,28,31,34,37,40,43,46,49
    };
    static constexpr double PHI_T = 51.2;
    static constexpr double PHI_R = 50.73;
    static constexpr double HTG = 10;
    static constexpr double HRG = 10;
    static constexpr double TX_GAIN = 20;
    static constexpr double RX_GAIN = 5;
    static constexpr double DN = 53;
    static constexpr double N0 = 328;
    static constexpr double DIST_COAST_TX = 500;
    static constexpr double DIST_COAST_RX = 500;
    static constexpr double DRY_PRESSURE_HPA = 1013;
    static constexpr double TEMP_C = 15;
    static constexpr ClutterModel::ClutterType TX_CLUTTER_TYPE = ClutterModel::ClutterType::IndustrialZone;
    static constexpr ClutterModel::ClutterType RX_CLUTTER_TYPE = ClutterModel::ClutterType::IndustrialZone;
    static constexpr auto POL = Enumerations::PolarizationType::VerticalPolarized;
    static constexpr double INPUT_LAT = (PHI_T+PHI_R)/2.0;
    static constexpr double TEMP_K = TEMP_C + 273.15;
};

#endif /* TEST_PROFILE_INDUSTRIAL_H */