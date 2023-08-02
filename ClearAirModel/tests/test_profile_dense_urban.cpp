#include "gtest/gtest.h"
#include "test_profile_dense_urban.h"
#include "ClearAirModel/ClearAirModelHelpers.h"
#include "ClearAirModel/BasicProp.h"
#include "ClearAirModel/DiffractionLoss.h"
#include "ClearAirModel/TropoScatter.h"
#include "ClearAirModel/AnomolousProp.h"
#include "ClearAirModel/P452TotalAttenuation.h"
#include <filesystem>

//Validation data from ITU validation spreadsheets in results folder
//in ITU validation dataset titled "Validation examples for software implementation of Recommendation ITU-R P.452"
//https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/R19-WP3M-C-0364!N18-P2!ZIP-E.zip

//Intermediate results from the transmission losses tab of the spreadsheets in the validation examples folder are inaccurate.
//They do not account for the 0.1dB difference in the free space path loss constant changed between 452-16 and 452-17. 
//As such, results from the spreadsheets in the results folder are preferred

//NOTE: this does not test the clutter model

namespace {
	// Use when expected an exact match
	double constexpr TOLERANCE_STRICT = 1.0e-6;
    // Use when different speed of light constant used
    // The validation data will match to a stricter tolerance if 2.998*1e8 is used for speed of light
    double constexpr TOLERANCE = 1.0e-3;
    const std::filesystem::path clearAirDataRelPath = "tests/test_paths";
	const std::filesystem::path clearAirDataFullPath = CMAKE_CLEARAIR_SRC_DIR / clearAirDataRelPath;
}

//allocate memory for the path
PathProfile::Path UrbanProfileTests::K_PATH;


TEST_F(UrbanProfileTests, calcP452TotalAttenuationTest){

    //basic loss prediction not exceeded for p percent
    const std::vector<double> EXPECTED_LOSS = {
        149.347206,100.3239608,100.5747582,102.5467334,107.8931615,
        123.9348437,140.3584782,144.4500525,147.9780989,151.5030596,
        155.0283871,158.5565199,162.0929341,165.6578413,169.5260939,
        173.0729193,176.9787614,179.1974968,147.6345642,148.2466785,
        148.4937753,148.6081691,148.6994773,148.7771238,148.8457268,
        148.9079392,148.9654408,149.0193735,149.0705578,149.1196107,
        149.1670157,149.2131659,149.2583925,149.302985,149.347206
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        const auto p452Model = ClearAirModel::p452_TotalAttenuation(
            FREQ_GHZ_LIST[freqInd],
            P_LIST[freqInd],
            K_PATH,
            HTG,
            HRG,
            INPUT_LAT,
            TX_GAIN,
            RX_GAIN,
            POL,
            DIST_COAST_TX,
            DIST_COAST_RX,
            DN,
            N0,
            TEMP_K,
            DRY_PRESSURE_HPA,
            TX_CLUTTER_TYPE,
            RX_CLUTTER_TYPE
        );
        double LOSS_VAL = p452Model.getTotalTransmissionLoss_dB();
        EXPECT_NEAR(EXPECTED_LOSS[freqInd],LOSS_VAL,TOLERANCE);
    }
}

