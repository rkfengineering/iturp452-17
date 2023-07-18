#include "gtest/gtest.h"
#include "EffectiveEarth.h"
#include <filesystem>

//Validation data from ITU validation spreadsheet titled "delB_valid_temp.xlsx", page "outputs"
//embedded in ITU validation document titled "Validation Examples for the delta Bullington diffraction prediction method"
//https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/Validation%20examples%20for%20the%20delta%20Bullington%20diffraction%20prediction%20method.docx

//If there's a way to use relative paths that would be nice. The paths from the validation data have been isolated as separate
//csv files and put in a subdirectory in the tests folder
const std::filesystem::path testPathFileDir("/home/ayeh/itu/ituModels/iturp452/ClearAirModel/tests/test_paths");

namespace {
	// Use when expected an exact match
	double constexpr TOLERANCE = 1.0e-6;
}

//From intermediate values in test_result_mixed_109.csv
TEST(EffectiveEarthTests, EffectiveEarthTests_effRadiusp50km){
    const double INPUT_DN = 53;
    const double EXPECTED_AE = 9617.759615;
	const double VAL_AE = EffectiveEarth::eff_radius_p50_km(INPUT_DN);

	EXPECT_NEAR(EXPECTED_AE,VAL_AE,TOLERANCE);
}

/*
TEST(BasicPropTests, BasicPropTests_pathLossWithGasAndMultipath){
	// Arrange
	const std::vector<double> FREQ_GHZ_LIST = {
		0.2,0.1,0.15,0.225,0.3375,0.50625,0.759375,1.1390625,
		1.70859375,2.562890625,3.844335938,5.766503906,8.649755859,
		12.97463379,19.46195068,29.19292603,43.78938904,50,0.2,0.2,
		0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2
	};
    const std::vector<double> HTS_LIST = {
        30,50,20,40,70,
        30,50,20,40,70,
        30,50,20,40,70,
        30,50,20,40,70
    };
    const std::vector<double> HRS_LIST = {
        30,10,20,50,5,
        30,10,20,50,5,
        30,10,20,50,5,
        30,10,20,50,5
    };
    const std::vector<double> FREQ_MHZ_LIST = {
        1000,2500,600,200,150,
        1000,2500,600,200,150,
        1000,2500,600,200,150,
        1000,2500,600,200,150
    };
    
    const std::vector<double> EXPECTED_HTSP = {
        211.3366061,233.4945304,204.7165548,215.727222,252.3138569,37.97366242,
        57.97664284,27.97366242,47.97366242,78.00115495,110.1894677,132.6734925,
        102.8893122,118.7202148,152.2011247,139.8495825,159.8495825,129.8495825,149.8495825,179.8495825
    };
    const std::vector<double> EXPECTED_HRSP = {
        30,10,20,50,5,119.2363492,99.82603022,109.2363492,139.2363492,99.67579634,30,
        11.9416848,22.52524312,50,5.664441408,543.1729418,523.1729418,533.1729418,
        563.1729418,518.1729418,
    };

    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
        //hts,hrs needs to be in m asl, frequency in GHz
        const PathProfile::Path p = PROFILE_LIST[PATH_LIST[pathInd]-1];
        const EffectiveEarth::TxRxPair EFF_HEIGHT = EffectiveEarth::smoothEarthHeights_diffractionModel(
            p,HTS_LIST[pathInd],HRS_LIST[pathInd]);
        const double HTS_AMSL = HTS_LIST[pathInd]+p.front().h_masl;
        const double HRS_AMSL = HRS_LIST[pathInd]+p.back().h_masl;
        //Equation 38 to calculate modified antenna heights
        EXPECT_NEAR(EXPECTED_HTSP[pathInd],HTS_AMSL-EFF_HEIGHT.tx_val,TOLERANCE);
        EXPECT_NEAR(EXPECTED_HRSP[pathInd],HRS_AMSL-EFF_HEIGHT.rx_val,TOLERANCE);
    }
}
*/