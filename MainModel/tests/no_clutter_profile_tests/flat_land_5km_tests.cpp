#include "gtest/gtest.h"
#include "flat_land_5km_testclass.h"
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
PathProfile::Path FlatLand5kmProfileTests::K_PATH;

//test height gain model heights
TEST_F(FlatLand5kmProfileTests, ClutterModelTests_calcHeightGainModelHeightsTest){
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
TEST_F(FlatLand5kmProfileTests, calcTimePercentBeta0Test){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EXPECTED_B0 = 7.005407788;
	EXPECT_NEAR(EXPECTED_B0,B0_PERCENT,TOLERANCE_STRICT);
}
TEST_F(FlatLand5kmProfileTests, calcFracOverSeaTest){
	// Arrange
    const double VAL_FRAC= K_PATH.calcFracOverSea();
    const double EXPECTED_FRAC = 0.0;
	EXPECT_NEAR(EXPECTED_FRAC,VAL_FRAC,TOLERANCE_STRICT);
}
TEST_F(FlatLand5kmProfileTests, calcLongestContiguousInlandDistanceTest){
	// Arrange
    const double VAL_DIST= K_PATH.calcLongestContiguousInlandDistance_km();
    const double EXPECTED_DIST = 5.0;
	EXPECT_NEAR(EXPECTED_DIST,VAL_DIST,TOLERANCE_STRICT);
}

TEST_F(FlatLand5kmProfileTests, BasicPropTests_calcPathLossWithGasAndMultipathTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LBFSG = {
        112.4345867,86.38041127,89.90342433,93.42765208,96.95382173,
        100.4821208,104.0112868,107.5393281,111.0654696,114.5904552,
        118.115811,121.6439947,125.1805265,128.7457812,132.6168275,
        136.1638546,140.0727936,142.3001315,112.4345867,112.4345867,
        112.4345867,112.4345867,112.4345867,112.4345867,112.4345867,
        112.4345867,112.4345867,112.4345867,112.4345867,112.4345867,
        112.4345867,112.4345867,112.4345867,112.4345867,112.4345867
    };
    const std::vector<double> EXPECTED_LB0P = {
        112.4256108,86.37143537,89.89444842,93.41867617,96.94484583,
        100.4731449,104.0023109,107.5303522,111.0564937,114.5814793,
        118.1068351,121.6350188,125.1715506,128.7368053,132.6078516,
        136.1548787,140.0638177,142.2911556,110.6965059,111.3124255,
        111.5610584,111.7195262,111.8360926,111.9283451,112.0046968,
        112.0698316,112.1266269,112.176978,112.2221992,112.2632399,
        112.3008081,112.3354458,112.3675773,112.3975409,112.4256108
    };
    const std::vector<double> EXPECTED_LB0B = {
        111.5614015,85.50722604,89.03023909,92.55446684,96.08063649,
        99.6089356,103.1381016,106.6661429,110.1922844,113.71727,
        117.2426258,120.7708094,124.3073412,127.872596,131.7436423,
        135.2906693,139.1996083,141.4269462,111.5614015,111.5614015,
        111.5614015,111.5614015,111.5614015,111.5614015,111.5614015,
        111.5614015,111.5614015,111.5614015,111.5614015,111.5614015,
        111.5614015,111.5614015,111.5614015,111.5614015,111.5614015
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

        ITUR_P452::CommonInputs commonInputs;
        commonInputs.freq_GHz=FREQ_GHZ_LIST[freqInd];
        commonInputs.p_percent = P_LIST[freqInd];
        commonInputs.height_tx_asl_m=height_tx_asl_m;
        commonInputs.height_rx_asl_m=height_rx_asl_m;
        commonInputs.path=mod_path;
        commonInputs.fracOverSea=SEA_FRAC;
        commonInputs.timePercentBeta0=B0_PERCENT;
        commonInputs.d_tot_km=mod_path.back().d_km;

        const auto BasicPropModel = ITUR_P452::BasicProp(commonInputs,TEMP_K, DRY_PRESSURE_HPA,HorizonDistances);

        double VAL_LBFSG,VAL_LB0P,VAL_LB0B;
	    BasicPropModel.calcTransmissionlosses_dB(VAL_LBFSG,VAL_LB0P,VAL_LB0B);

        //WARNING a higher tolerance is used. 
        EXPECT_NEAR(EXPECTED_LBFSG[freqInd],VAL_LBFSG,TOLERANCE);
        EXPECT_NEAR(EXPECTED_LB0P[freqInd],VAL_LB0P,TOLERANCE);
        EXPECT_NEAR(EXPECTED_LB0B[freqInd],VAL_LB0B,TOLERANCE);
    }
}

TEST_F(FlatLand5kmProfileTests, DiffractionLossTests_calcDiffractionLossTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LD50 = {
        0,19.99720128,16.69932053,13.48785052,10.32414568,7.16721976,
        3.9608281,0.59152715,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    const std::vector<double> EXPECTED_LDP = {
        0,19.99493392,16.69677175,13.48498097,10.32090779,7.1635521,
        3.95664252,0.58673311,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        //apply height gain clutter model
        const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ_LIST[freqInd], K_PATH, 
                HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
        const auto mod_path = ClutterResults.modifiedPath;
        const auto [hg_height_tx_m, hg_height_rx_m] = ClutterResults.modifiedHeights_m;

        const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
        const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

        ITUR_P452::CommonInputs commonInputs;
        commonInputs.freq_GHz=FREQ_GHZ_LIST[freqInd];
        commonInputs.p_percent = P_LIST[freqInd];
        commonInputs.height_tx_asl_m=height_tx_asl_m;
        commonInputs.height_rx_asl_m=height_rx_asl_m;
        commonInputs.path=mod_path;
        commonInputs.fracOverSea=SEA_FRAC;
        commonInputs.timePercentBeta0=B0_PERCENT;
        commonInputs.d_tot_km=mod_path.back().d_km;

        const auto DiffractionModel = ITUR_P452::DiffractionLoss(
            commonInputs,
            DN,
            POL
        );

        double diffractionLoss_median_dB, diffractionLoss_p_percent_dB;
	    DiffractionModel.calcDiffractionLoss_dB(diffractionLoss_median_dB, diffractionLoss_p_percent_dB);
        EXPECT_NEAR(EXPECTED_LD50[freqInd],diffractionLoss_median_dB,TOLERANCE);
        EXPECT_NEAR(EXPECTED_LDP[freqInd],diffractionLoss_p_percent_dB,TOLERANCE);
    }
}

TEST_F(FlatLand5kmProfileTests, TroposcatterLossTests_calcTroposcatterLossTest){
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
        162.1737211,125.3835186,130.8549714,136.1725987,141.3371192,
        146.3486964,151.2060373,155.9071302,160.4511464,164.8386521,
        169.0707337,173.1487884,177.0754523,180.8596085,184.6152725,
        188.0997816,191.7337422,193.7777318,147.9036526,151.7650386,
        153.4979669,154.6802706,155.600871,156.3682457,157.0356905,
        157.633772,158.182017,158.6940484,159.1802158,159.6491798,
        160.1091358,160.5692078,161.0419885,161.5508925,162.1737211
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        //input effective antanna heights relative to ground, frequency in GHz      

        ITUR_P452::CommonInputs commonInputs;
        commonInputs.freq_GHz=FREQ_GHZ_LIST[freqInd];
        commonInputs.p_percent = P_LIST[freqInd];
        commonInputs.height_tx_asl_m=height_tx_asl_m;
        commonInputs.height_rx_asl_m=height_rx_asl_m;
        commonInputs.path=mod_path;//not used
        commonInputs.fracOverSea=0;//not used
        commonInputs.timePercentBeta0=3;//not used
        commonInputs.d_tot_km=D_TOT_KM;

        const double LOSS_VAL = ITUR_P452::TropoScatter::calcTroposcatterLoss_dB(
            commonInputs,
            HorizonAngles,
            EFF_RADIUS_MED_KM,
            N0,
            TX_GAIN,
            RX_GAIN,
            TEMP_K,
            DRY_PRESSURE_HPA
        );
        EXPECT_NEAR(EXPECTED_LBS[freqInd],LOSS_VAL,TOLERANCE);
    }
}
//check horizon elevation angle and distance calculations using mixed terrain path
//line of sight path tested
TEST_F(FlatLand5kmProfileTests, calcHorizonAnglesAndDistances_LineOfSightTest){
    const double FREQ_GHZ = 2.0;
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);

    const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
        const auto mod_path = ClutterResults.modifiedPath;
    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterResults.modifiedHeights_m;

    const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
    const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

    const double EXPECTED_THETA_T = -0.25993579;
    const double EXPECTED_THETA_R = -0.25993579;
    const double EXPECTED_DLT = 2.5;
    const double EXPECTED_DLR = 2.5;

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
TEST_F(FlatLand5kmProfileTests, AnomalousProp_calcAnomalousPropLossTest){

    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    //basic transmission loss from ducting/layer refraction
    const std::vector<double> EXPECTED_LBA = {
        184.1412539,190.6870785,188.5163415,184.3671318,178.3343171,
        172.1887881,175.717954,179.2459954,182.7721369,186.2971224,
        189.8224782,193.3506619,196.8871937,200.4524484,204.3234947,
        207.8705218,211.7794608,214.0067987,111.4343185,117.3315672,
        122.4756036,127.3549699,132.0792493,136.695792,141.2305881,
        145.6999065,150.1148228,154.4833313,158.8114568,163.1038922,
        167.3643906,171.5960177,175.8013223,179.9824551,184.1412539
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

        ITUR_P452::CommonInputs commonInputs;
        commonInputs.freq_GHz=FREQ_GHZ_LIST[freqInd];
        commonInputs.p_percent = P_LIST[freqInd];
        commonInputs.height_tx_asl_m=height_tx_asl_m;
        commonInputs.height_rx_asl_m=height_rx_asl_m;
        commonInputs.path=mod_path;
        commonInputs.fracOverSea=SEA_FRAC;
        commonInputs.timePercentBeta0=B0_PERCENT;
        commonInputs.d_tot_km=mod_path.back().d_km;

        const auto AnomalousModel = ITUR_P452::AnomalousProp(
            commonInputs,
            TEMP_K,
            DRY_PRESSURE_HPA,
            DIST_COAST_TX,
            DIST_COAST_RX,
            EFF_RADIUS_MED_KM,
            HORIZON_VALS
        );
        const double LOSS_VAL = AnomalousModel.calcAnomalousPropLoss_dB();
        EXPECT_NEAR(EXPECTED_LBA[freqInd],LOSS_VAL,TOLERANCE);
    }
}

TEST_F(FlatLand5kmProfileTests, calcP452TotalAttenuationTest){
    //basic loss prediction not exceeded for p percent
    const std::vector<double> EXPECTED_LOSS = {
        112.4197947,106.3624412,106.5878792,106.9006589,107.2631202,
        107.6344864,107.9572519,108.115982,111.0506776,114.5756632,
        118.101019,121.6292026,125.1657344,128.7309891,132.6020354,
        136.1490625,140.0580015,142.2853394,110.6965059,111.3124255,
        111.5610584,111.6761634,111.7680392,111.8461684,111.9151979,
        111.977797,112.0356561,112.0899242,112.1414266,112.1907845,
        112.2384842,112.2849212,112.330429,112.3752988,112.4197947
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

