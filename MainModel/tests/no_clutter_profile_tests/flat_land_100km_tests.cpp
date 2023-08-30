#include "gtest/gtest.h"
#include "flat_land_100km_testclass.h"
#include "MainModel/ClearAirModelHelpers.h"
#include "MainModel/BasicProp.h"
#include "MainModel/DiffractionLoss.h"
#include "MainModel/TropoScatter.h"
#include "MainModel/AnomalousProp.h"
#include "MainModel/P452TotalAttenuation.h"
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
PathProfile::Path FlatLand100kmProfileTests::K_PATH;

//test height gain model heights
TEST_F(FlatLand100kmProfileTests, ClutterModelTests_calcHeightGainModelHeightsTest){
    const double FREQ_GHZ = 2;
    const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
    const auto HEIGHT_GAIN_RESULTS = ClutterResults.modifiedHeights_m;
    const auto MOD_PATH = ClutterResults.modifiedPath;
    const double EXPECTED_HTS = 10.0;
    const double EXPECTED_HRS = 10.0;
    EXPECT_NEAR(EXPECTED_HTS,HEIGHT_GAIN_RESULTS.first + MOD_PATH.front().h_asl_m,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_HRS,HEIGHT_GAIN_RESULTS.second+ MOD_PATH.back().h_asl_m,TOLERANCE_STRICT);
}
TEST_F(FlatLand100kmProfileTests, calcTimePercentBeta0Test){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EXPECTED_B0 = 1.224161171;
	EXPECT_NEAR(EXPECTED_B0,B0_PERCENT,TOLERANCE_STRICT);
}
TEST_F(FlatLand100kmProfileTests, calcFracOverSeaTest){
	// Arrange
    const double VAL_FRAC= K_PATH.calcFracOverSea();
    const double EXPECTED_FRAC = 0.0;
	EXPECT_NEAR(EXPECTED_FRAC,VAL_FRAC,TOLERANCE_STRICT);
}
TEST_F(FlatLand100kmProfileTests, calcLongestContiguousInlandDistanceTest){
	// Arrange
    const double VAL_DIST= K_PATH.calcLongestContiguousInlandDistance_km();
    const double EXPECTED_DIST = 100.0;
	EXPECT_NEAR(EXPECTED_DIST,VAL_DIST,TOLERANCE_STRICT);
}

TEST_F(FlatLand100kmProfileTests, BasicPropTests_calcPathLossWithGasAndMultipathTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LBFSG = {
        139.1123342,112.4202238,115.9658064,119.5356829,123.1443976,
        126.7957013,130.4643423,134.1104904,137.7186421,141.3036749,
        144.8961128,148.5451073,152.3610648,156.7514809,167.2577284,
        171.2835915,182.5476933,205.2060261,139.1123342,139.1123342,
        139.1123342,139.1123342,139.1123342,139.1123342,139.1123342,
        139.1123342,139.1123342,139.1123342,139.1123342,139.1123342,
        139.1123342,139.1123342,139.1123342,139.1123342,139.1123342
    };
    const std::vector<double> EXPECTED_LB0P = {
        139.0909092,112.3987988,115.9443814,119.5142579,123.1229726,
        126.7742763,130.4429173,134.0890654,137.6972171,141.2822499,
        144.8746878,148.5236824,152.3396399,156.7300559,167.2363034,
        171.2621665,182.5262683,205.1846011,134.9636298,136.4337964,
        137.0272695,137.4055236,137.6837612,137.903963,138.0862103,
        138.2416835,138.3772509,138.4974361,138.6053767,138.7033387,
        138.7930119,138.8756902,138.9523864,139.0239079,139.0909092
    };
    const std::vector<double> EXPECTED_LB0B = {
        135.1781224,108.486012,112.0315946,115.6014711,119.2101858,
        122.8614895,126.5301305,130.1762786,133.7844303,137.3694632,
        140.961901,144.6108956,148.4268531,152.8172692,163.3235166,
        167.3493797,178.6134815,201.2718144,135.1781224,135.1781224,
        135.1781224,135.1781224,135.1781224,135.1781224,135.1781224,
        135.1781224,135.1781224,135.1781224,135.1781224,135.1781224,
        135.1781224,135.1781224,135.1781224,135.1781224,135.178122
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        //apply height gain clutter model
        const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ_LIST[freqInd], K_PATH, 
                HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
        const auto mod_path = ClutterResults.modifiedPath;
        const auto [hg_height_tx_m, hg_height_rx_m] = ClutterResults.modifiedHeights_m;

        const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
        const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

        const auto [HorizonAngles, HorizonDistances] = ITUR_P452::Helpers::calcHorizonAnglesAndDistances(
            mod_path, height_tx_asl_m, height_rx_asl_m, EFF_RADIUS_MED_KM, FREQ_GHZ_LIST[freqInd]
        );

        const auto BasicPropModel = ITUR_P452::BasicProp(mod_path.back().d_km, height_tx_asl_m, height_rx_asl_m, 
                FREQ_GHZ_LIST[freqInd],TEMP_K, DRY_PRESSURE_HPA, SEA_FRAC,P_LIST[freqInd],B0_PERCENT,HorizonDistances);

        double VAL_LBFSG,VAL_LB0P,VAL_LB0B;
	    BasicPropModel.calcTransmissionlosses_dB(VAL_LBFSG,VAL_LB0P,VAL_LB0B);

        //WARNING a higher tolerance is used. 
        EXPECT_NEAR(EXPECTED_LBFSG[freqInd],VAL_LBFSG,TOLERANCE);
        EXPECT_NEAR(EXPECTED_LB0P[freqInd],VAL_LB0P,TOLERANCE);
        EXPECT_NEAR(EXPECTED_LB0B[freqInd],VAL_LB0B,TOLERANCE);
    }
}

TEST_F(FlatLand100kmProfileTests, DiffractionLossTests_calcDiffractionLossTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LD50 = {
        93.39174946,65.52126845,65.9425145,67.18446656,69.36054332,
        72.59738582,77.03445261,82.82216141,90.11777129,99.07834171,
        109.8513632,122.5674981,137.3473332,154.4973108,174.5237975,
        197.797126,224.7109561,234.3818887,93.39174946,93.39174946,
        93.39174946,93.39174946,93.39174946,93.39174946,93.39174946,
        93.39174946,93.39174946,93.39174946,93.39174946,93.39174946,
        93.39174946,93.39174946,93.39174946,93.39174946,93.39174946
    };
    const std::vector<double> EXPECTED_LDP = {
        93.02369321,65.42513337,65.82434308,67.03965383,69.1852776,
        72.38733792,76.78472275,82.52723215,89.77147296,98.67383494,
        109.3811299,122.0232939,136.7199801,153.7744116,173.6901643,
        196.8372323,223.6072297,233.2268493,60.26607136,67.61036241,
        71.65920745,74.52064635,76.80670136,78.75071337,80.46830742,
        82.02590048,83.4655513,84.8158494,86.09733467,87.32545838,
        88.51232326,89.66777021,90.8000948,91.91654587,93.02369321
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        //apply height gain clutter model
        const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ_LIST[freqInd], K_PATH, 
                HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
        const auto mod_path = ClutterResults.modifiedPath;
        const auto [hg_height_tx_m, hg_height_rx_m] = ClutterResults.modifiedHeights_m;

        const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
        const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

        const auto DiffractionModel = ITUR_P452::DiffractionLoss(
            mod_path, 
            height_tx_asl_m,
            height_rx_asl_m,
            FREQ_GHZ_LIST[freqInd],
            DN,
            POL,
            P_LIST[freqInd],
            B0_PERCENT,
            SEA_FRAC
        );

        double diffractionLoss_median_dB, diffractionLoss_p_percent_dB;
	    DiffractionModel.calcDiffractionLoss_dB(diffractionLoss_median_dB, diffractionLoss_p_percent_dB);
        EXPECT_NEAR(EXPECTED_LD50[freqInd],diffractionLoss_median_dB,TOLERANCE);
        EXPECT_NEAR(EXPECTED_LDP[freqInd],diffractionLoss_p_percent_dB,TOLERANCE);
    }
}

TEST_F(FlatLand100kmProfileTests, TroposcatterLossTests_calcTroposcatterLossTest){
	// Arrange
    const double FREQ_GHZ = 2.0;
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);
    //apply height gain clutter model
    const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
    const auto mod_path = ClutterResults.modifiedPath;
    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterResults.modifiedHeights_m;

    const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
    const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;
    const double D_TOT_KM = mod_path.back().d_km;

    const auto [HorizonAngles, HorizonDistances] = ITUR_P452::Helpers::calcHorizonAnglesAndDistances(
        mod_path,
        height_tx_asl_m,
        height_rx_asl_m,
        EFF_RADIUS_MED_KM,
        FREQ_GHZ
    );

    //basic transmission loss from troposcatter
    const std::vector<double> EXPECTED_LBS = {
        193.1410293,155.7284828,161.2225754,166.5859317,171.832924,
        176.9668217,181.9617666,186.7775244,191.3975216,195.8330798,
        200.1059295,204.244014,208.3000533,212.4517121,218.9792996,
        223.02956,233.0146254,254.025963,178.8709609,182.7323469,
        184.4652752,185.6475788,186.5681792,187.3355539,188.0029987,
        188.6010803,189.1493252,189.6613566,190.147524,190.6164881,
        191.076444,191.536516,192.0092967,192.5182008,193.1410293
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        //input effective antanna heights relative to ground, frequency in GHz        
        const double LOSS_VAL = ITUR_P452::TropoScatter::calcTroposcatterLoss_dB(
            D_TOT_KM, 
            FREQ_GHZ_LIST[freqInd],
            height_tx_asl_m,
            height_rx_asl_m,
            HorizonAngles,
            EFF_RADIUS_MED_KM,
            N0,
            TX_GAIN,
            RX_GAIN,
            TEMP_K,
            DRY_PRESSURE_HPA,
            P_LIST[freqInd]
        );
        EXPECT_NEAR(EXPECTED_LBS[freqInd],LOSS_VAL,TOLERANCE);
    }
}
//check horizon elevation angle and distance calculations using mixed terrain path
//line of sight path tested
TEST_F(FlatLand100kmProfileTests, calcHorizonAnglesAndDistances_LineOfSightTest){
    const double FREQ_GHZ = 2.0;
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);

    const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
        const auto mod_path = ClutterResults.modifiedPath;
    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterResults.modifiedHeights_m;

    const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
    const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

    const double EXPECTED_THETA_T = -1.442104943;
    const double EXPECTED_THETA_R = -1.442104943;
    const double EXPECTED_DLT = 14.0;
    const double EXPECTED_DLR = 14.0;

    const auto [HorizonAngles, HorizonDistances] = ITUR_P452::Helpers::calcHorizonAnglesAndDistances(
        mod_path,
        height_tx_asl_m,
        height_rx_asl_m,
        EFF_RADIUS_MED_KM,
        FREQ_GHZ
    );
    const auto [horizonElevation_tx_mrad, horizonElevation_rx_mrad] = HorizonAngles;
    const auto [horizonDist_tx_km, horizonDist_rx_km] = HorizonDistances;

    EXPECT_NEAR(EXPECTED_THETA_T,horizonElevation_tx_mrad,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_THETA_R,horizonElevation_rx_mrad,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_DLT,horizonDist_tx_km,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_DLR,horizonDist_rx_km,TOLERANCE_STRICT);
}
//no direct unit test for the helpers yet. they're all lumped together in this validation data check
TEST_F(FlatLand100kmProfileTests, AnomalousProp_calcAnomalousPropLossTest){

    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    //basic transmission loss from ducting/layer refraction
    const std::vector<double> EXPECTED_LBA = {
        236.8785224,239.9113088,238.0058293,234.1800766,228.547818,
        222.8893253,226.9746788,231.0978435,235.2520429,239.4621444,
        243.7701072,248.2381733,252.9917337,258.4554372,270.1902921,
        275.6225595,288.4965925,311.7305693,152.4407703,169.3239336,
        179.180923,186.6998711,192.9512552,198.3831106,203.2319652,
        207.6401305,211.7008439,215.478985,219.0217798,222.3648269,
        225.5357215,228.5563501,231.4444063,234.2144269,236.8785224
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ_LIST[freqInd], K_PATH, 
                HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
        const auto mod_path = ClutterResults.modifiedPath;
        const auto [hg_height_tx_m, hg_height_rx_m] = ClutterResults.modifiedHeights_m;

        const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
        const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

        const auto HORIZON_VALS = ITUR_P452::Helpers::calcHorizonAnglesAndDistances(
            mod_path, height_tx_asl_m, height_rx_asl_m, EFF_RADIUS_MED_KM, FREQ_GHZ_LIST[freqInd]);
        const auto AnomalousModel = ITUR_P452::AnomalousProp(
            mod_path,
            FREQ_GHZ_LIST[freqInd],
            height_tx_asl_m,
            height_rx_asl_m,
            TEMP_K,
            DRY_PRESSURE_HPA,
            DIST_COAST_TX,
            DIST_COAST_RX,
            P_LIST[freqInd],
            B0_PERCENT,
            EFF_RADIUS_MED_KM,
            HORIZON_VALS,
            SEA_FRAC
        );
        const double LOSS_VAL = AnomalousModel.calcAnomalousPropLoss_dB();
        EXPECT_NEAR(EXPECTED_LBA[freqInd],LOSS_VAL,TOLERANCE);
    }
}

TEST_F(FlatLand100kmProfileTests, calcP452TotalAttenuationTest){
    //basic loss prediction not exceeded for p percent
    const std::vector<double> EXPECTED_LOSS = {
        193.1410293,155.7284001,161.2224066,166.5857114,171.8327495,
        176.9667427,181.9617474,186.7775221,191.3975214,195.8330798,
        200.1059295,204.244014,208.3000533,212.4517121,218.9792996,
        223.02956,233.0146254,254.025963,152.4433215,169.3196367,
        178.9984917,184.6055058,186.456268,187.3221915,188.0010455,
        188.6007423,189.1492581,189.6613417,190.1475204,190.6164871,
        191.0764437,191.5365159,192.0092967,192.5182007,193.1410293
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        const auto p452Model = ITUR_P452::TotalClearAirAttenuation(
            FREQ_GHZ_LIST[freqInd],
            P_LIST[freqInd],
            K_PATH,
            HTG,
            HRG,
            INPUT_LAT,
            DN,
            TX_CLUTTER_TYPE,
            RX_CLUTTER_TYPE
        );
        const double LOSS_VAL = p452Model.calcTotalClearAirAttenuation(
            TEMP_K,
            DRY_PRESSURE_HPA,
            DIST_COAST_TX,
            DIST_COAST_RX,
            N0,
            TX_GAIN,
            RX_GAIN,
            POL
        );
        EXPECT_NEAR(EXPECTED_LOSS[freqInd],LOSS_VAL,TOLERANCE);
    }
}

