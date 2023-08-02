#include "gtest/gtest.h"
#include "test_profile_dense_suburban.h"
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
PathProfile::Path SuburbanProfileTests::K_PATH;

//test height gain model heights
TEST_F(SuburbanProfileTests, ClutterLossTests_calcHeightGainModelHeightsTest){
    const double FREQ_GHZ = 2;
    const auto ClutterModel = ClearAirModel::ClutterLoss(FREQ_GHZ, K_PATH, 
        HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);

    const auto HEIGHT_GAIN_RESULTS = ClutterModel.getHeightGainModelHeights_m();
    const auto MOD_PATH = ClutterModel.getModifiedPath();

    const double EXPECTED_HTS = 12.0;//above sea level
    const double EXPECTED_HRS = 12.0;
    EXPECT_NEAR(EXPECTED_HTS,HEIGHT_GAIN_RESULTS.first + MOD_PATH.front().h_asl_m,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_HRS,HEIGHT_GAIN_RESULTS.second+ MOD_PATH.back().h_asl_m,TOLERANCE_STRICT);
}
TEST_F(SuburbanProfileTests, calcTimePercentBeta0Test){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EXPECTED_B0 = 7.005407788;
	EXPECT_NEAR(EXPECTED_B0,B0_PERCENT,TOLERANCE_STRICT);
}
TEST_F(SuburbanProfileTests, calcFracOverSeaTest){
	// Arrange
    const double VAL_FRAC= K_PATH.calcFracOverSea();
    const double EXPECTED_FRAC = 0.0;
	EXPECT_NEAR(EXPECTED_FRAC,VAL_FRAC,TOLERANCE_STRICT);
}
TEST_F(SuburbanProfileTests, calcLongestContiguousInlandDistanceTest){
	// Arrange
    const double VAL_DIST= K_PATH.calcLongestContiguousInlandDistance_km();
    const double EXPECTED_DIST = 5.0;
	EXPECT_NEAR(EXPECTED_DIST,VAL_DIST,TOLERANCE_STRICT);
}

TEST_F(SuburbanProfileTests, BasicPropTests_calcPathLossWithGasAndMultipathTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ClearAirModel::ClearAirModelHelpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LBFSG = {
        112.3645435,86.31063663,89.83364018,93.35784871,96.88398361,
        100.4122309,103.9413382,107.4693298,110.9954367,114.520397,
        118.0457246,121.5738574,125.1102715,128.6751788,132.5434313,
        136.0902568,139.9960989,142.2148343,112.3645435,112.3645435,
        112.3645435,112.3645435,112.3645435,112.3645435,112.3645435,
        112.3645435,112.3645435,112.3645435,112.3645435,112.3645435,
        112.3645435,112.3645435,112.3645435,112.3645435,112.3645435
    };
    const std::vector<double> EXPECTED_LB0P = {
        112.355623,86.30171618,89.82471973,93.34892826,96.87506316,
        100.4033105,103.9324177,107.4604093,110.9865163,114.5114766,
        118.0368042,121.5649369,125.1013511,128.6662584,132.5345109,
        136.0813364,139.9871785,142.2059138,110.6372011,111.2493154,
        111.4964121,111.6539009,111.769747,111.8614296,111.9373096,
        112.0020419,112.0584864,112.1085263,112.1534682,112.1942553,
        112.2315914,112.2660151,112.2979481,112.3277266,112.355623
    };
    const std::vector<double> EXPECTED_LB0B = {
        111.4967531,85.44284621,88.96584976,92.49005829,96.01619319,
        99.5444405,103.0735478,106.6015394,110.1276463,113.6526066,
        117.1779342,120.706067,124.2424811,127.8073884,131.6756409,
        135.2224664,139.1283085,141.3470439,111.4967531,111.4967531,
        111.4967531,111.4967531,111.4967531,111.4967531,111.4967531,
        111.4967531,111.4967531,111.4967531,111.4967531,111.4967531,
        111.4967531,111.4967531,111.4967531,111.4967531,111.4967531
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        //apply height gain clutter model
        const auto ClutterModel = ClearAirModel::ClutterLoss(FREQ_GHZ_LIST[freqInd], K_PATH, 
                HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
        const auto mod_path = ClutterModel.getModifiedPath();
        const auto [hg_height_tx_m, hg_height_rx_m] = ClutterModel.getHeightGainModelHeights_m();
        const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
        const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

        const auto [HorizonAngles, HorizonDistances] = ClearAirModel::ClearAirModelHelpers::calcHorizonAnglesAndDistances(
            mod_path, height_tx_asl_m, height_rx_asl_m, EFF_RADIUS_MED_KM, FREQ_GHZ_LIST[freqInd]
        );

        const auto BasicPropModel = ClearAirModel::BasicProp(mod_path.back().d_km, height_tx_asl_m, height_rx_asl_m, 
                FREQ_GHZ_LIST[freqInd],TEMP_K, DRY_PRESSURE_HPA, SEA_FRAC,P_LIST[freqInd],B0_PERCENT,HorizonDistances);

        const double VAL_LBFSG = BasicPropModel.getFreeSpaceWithGasLoss_dB();
        const double VAL_LB0P = BasicPropModel.getBasicTransmissionLoss_p_percent_dB();
        const double VAL_LB0B = BasicPropModel.getBasicTransmissionLoss_b0_percent_dB();

        EXPECT_NEAR(EXPECTED_LBFSG[freqInd],VAL_LBFSG,TOLERANCE_STRICT);
        EXPECT_NEAR(EXPECTED_LB0P[freqInd],VAL_LB0P,TOLERANCE_STRICT);
        EXPECT_NEAR(EXPECTED_LB0B[freqInd],VAL_LB0B,TOLERANCE_STRICT);
    }
}

TEST_F(SuburbanProfileTests, DiffractionLossTests_calcDiffractionLossTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LD50 = {
        0,16.94131569,13.70404793,10.52393013,7.35075448,
        4.13851198,0.77100716,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    const std::vector<double> EXPECTED_LDP = {
        0,16.93924131,13.70171403,10.52129766,7.34777559,
        4.13511496,0.76711705,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        //apply height gain clutter model
        const auto ClutterModel = ClearAirModel::ClutterLoss(FREQ_GHZ_LIST[freqInd], K_PATH, 
                HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
        const auto mod_path = ClutterModel.getModifiedPath();
        const auto [hg_height_tx_m, hg_height_rx_m] = ClutterModel.getHeightGainModelHeights_m();
        const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
        const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

        const auto DiffractionModel = ClearAirModel::DiffractionLoss(
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

        const double diffractionLoss_median_dB = DiffractionModel.getDiffractionLoss_median_dB();
        const double diffractionLoss_p_percent_dB = DiffractionModel.getDiffractionLoss_p_percent_dB();

        EXPECT_NEAR(EXPECTED_LD50[freqInd],diffractionLoss_median_dB,TOLERANCE);
        EXPECT_NEAR(EXPECTED_LDP[freqInd],diffractionLoss_p_percent_dB,TOLERANCE);
    }
}

TEST_F(SuburbanProfileTests, TroposcatterLossTests_calcTroposcatterLossTest){
	// Arrange
    const double FREQ_GHZ = 2.0;
    const double EFF_RADIUS_MED_KM = ClearAirModel::ClearAirModelHelpers::calcMedianEffectiveRadius_km(DN);
    //apply height gain clutter model
    const auto ClutterModel = ClearAirModel::ClutterLoss(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
    const auto mod_path = ClutterModel.getModifiedPath();
    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterModel.getHeightGainModelHeights_m();
    const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
    const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;
    const double D_TOT_KM = mod_path.back().d_km;

    const auto [HorizonAngles, HorizonDistances] = ClearAirModel::ClearAirModelHelpers::calcHorizonAnglesAndDistances(
        mod_path,
        height_tx_asl_m,
        height_rx_asl_m,
        EFF_RADIUS_MED_KM,
        FREQ_GHZ
    );

    //basic transmission loss from troposcatter
    const std::vector<double> EXPECTED_LBS = {
        162.1036844,125.3137439,130.7851872,136.1027953,141.267281,
        146.2788067,151.1360896,155.8371343,160.3811185,164.768604,
        169.0006684,173.0786979,177.0053073,180.7893087,184.5438056,
        188.0280764,191.659363,193.6953663,147.8336159,151.6950019,
        153.4279302,154.6102339,155.5308342,156.298209,156.9656538,
        157.5637353,158.1119803,158.6240117,159.1101791,159.5791431,
        160.0390991,160.4991711,160.9719518,161.4808558,162.1036844
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        //input effective antanna heights relative to ground, frequency in GHz        
        const double LOSS_VAL = ClearAirModel::TropoScatter::calcTroposcatterLoss_dB(
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
TEST_F(SuburbanProfileTests, calcHorizonAnglesAndDistances_LineOfSightTest){
    const double FREQ_GHZ = 2.0;
    const double EFF_RADIUS_MED_KM = ClearAirModel::ClearAirModelHelpers::calcMedianEffectiveRadius_km(DN);

    const auto ClutterModel = ClearAirModel::ClutterLoss(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
    const auto mod_path = ClutterModel.getModifiedPath();
    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterModel.getHeightGainModelHeights_m();
    const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
    const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

    const double EXPECTED_THETA_T = -0.257856304;
    const double EXPECTED_THETA_R = -0.257856304;
    const double EXPECTED_DLT = 2.48;
    const double EXPECTED_DLR = 2.48;

    const auto [HorizonAngles, HorizonDistances] = ClearAirModel::ClearAirModelHelpers::calcHorizonAnglesAndDistances(
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
TEST_F(SuburbanProfileTests, AnomolousProp_calcAnomolousPropLossTest){

    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ClearAirModel::ClearAirModelHelpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    //basic transmission loss from ducting/layer refraction
    const std::vector<double> EXPECTED_LBA = {
        184.0711264,190.6172195,188.4464731,184.2972441,178.2643947,
        172.1188138,175.6479211,179.1759127,182.7020197,186.2269799,
        189.7523075,193.2804403,196.8168544,200.3817617,204.2500143,
        207.7968397,211.7026818,213.9214172,111.3643993,117.2615588,
        122.4055604,127.2849055,132.0091701,136.6257018,141.1604894,
        145.6298011,150.0447121,154.4132163,158.7413383,163.033771,
        167.2942673,171.5258927,175.7311961,179.9123281,184.0711264
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        const auto ClutterModel = ClearAirModel::ClutterLoss(FREQ_GHZ_LIST[freqInd], K_PATH, 
                HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
        const auto mod_path = ClutterModel.getModifiedPath();
        const auto [hg_height_tx_m, hg_height_rx_m] = ClutterModel.getHeightGainModelHeights_m();
        const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
        const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

        const auto HORIZON_VALS = ClearAirModel::ClearAirModelHelpers::calcHorizonAnglesAndDistances(
            mod_path, height_tx_asl_m, height_rx_asl_m, EFF_RADIUS_MED_KM, FREQ_GHZ_LIST[freqInd]);
        const auto AnomolousModel = ClearAirModel::AnomolousProp(
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
        const double LOSS_VAL = AnomolousModel.getAnomolousPropLoss_dB();
        EXPECT_NEAR(EXPECTED_LBA[freqInd],LOSS_VAL,TOLERANCE);
    }
}

TEST_F(SuburbanProfileTests, calcP452TotalAttenuationTest){

    //basic loss prediction not exceeded for p percent
    const std::vector<double> EXPECTED_LOSS = {
        114.7384424,103.3449365,103.6370257,104.0055519,104.5058701,
        105.8349089,107.0403958,109.8430716,113.3693356,116.8942959,
        120.4196235,123.9477563,127.4841704,131.0490777,134.9173303,
        138.4641557,142.3699978,144.5887332,113.0258006,113.6379149,
        113.8850116,113.9994055,114.0907137,114.1683601,114.2369632,
        114.2991756,114.3566772,114.4106099,114.4617942,114.5108471,
        114.5582521,114.6044023,114.6496288,114.6942214,114.7384424
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

