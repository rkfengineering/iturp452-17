
#ifndef DELTA_BULLINGTON_TESTS_H
#define DELTA_BULLINGTON_TESTS_H

#include "gtest/gtest.h"

#include <ClearAirModel/PathProfile.h>

#include <filesystem>

class DeltaBullingtonTests : public testing::Test {
protected:
	DeltaBullingtonTests() {}

	// Per-test-suite set-up.
	// Called before the first test in this test suite.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
        std::filesystem::path clearAirDataFullPath = CMAKE_CLEARAIR_SRC_DIR;
        clearAirDataFullPath /= std::filesystem::path("tests/test_paths");

        m_profile_list = {
            PathProfile::Path(clearAirDataFullPath/std::filesystem::path("dbull_path1.csv")),
            PathProfile::Path(clearAirDataFullPath/std::filesystem::path("dbull_path2.csv")),
            PathProfile::Path(clearAirDataFullPath/std::filesystem::path("dbull_path3.csv")),
            PathProfile::Path(clearAirDataFullPath/std::filesystem::path("dbull_path4.csv")),
        };
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
    static std::vector<PathProfile::Path> m_profile_list;

    //Assumptions from "Notes" tab of spreadsheet
    const double m_effEarthRadius_km = 8500.0;
    //Validation data uses 3e8 for speed of light
    //this implementation uses a slightly different constant
    const double m_deltaBullingtonTolerance_dB = 0.1;

};

#endif /* DELTA_BULLINGTON_TESTS_H */