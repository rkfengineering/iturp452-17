#include "gtest/gtest.h"

#include "P452/P452.h"

//Test to make sure interface works without runtime errors
//check to see if its greater than FSPL

namespace {
	// Use when expected an exact match
	double constexpr TOLERANCE = 1.0e-3;
}

TEST(P452WrapperTests, rawElevationInputTest){

	const std::vector<double> ELEVATION_LIST_M = {
		62.0, 62.0, 60.0, 66.0, 73.0, 88.0, 96.0, 108.0, 105.0, 84.0, 
        78.0, 63.0, 34.0, 38.0, 27.0, 19.0, 1.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	};
    const double stepDistance_km = 0.994291;
    const double txHeight_m = 10.0;
    const double rxHeight_m = 10.0;

    const double midpoint_lat_deg = 29.0002;
    const double midpoint_lon_deg = 48.25;

    const double freq_GHz = 0.3;
    const double timePercent = 50.0;

    const double EXPECTED_LOSS = 146.409;
    const double FSPL = 115.747;

    const double RES_LOSS = P452::calculateP452Loss_dB(txHeight_m, rxHeight_m, ELEVATION_LIST_M, stepDistance_km, 
            midpoint_lat_deg, midpoint_lon_deg, freq_GHz, timePercent);

    EXPECT_NEAR(EXPECTED_LOSS, RES_LOSS, TOLERANCE);
    EXPECT_TRUE(RES_LOSS>FSPL);
}