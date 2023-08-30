#include "gtest/gtest.h"
#include "dense_urban_testclass.h"
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

//test height gain model heights
TEST_F(UrbanProfileTests, ClutterModelTests_calcHeightGainModelHeightsTest){
    const double FREQ_GHZ = 2;
    const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
    const auto HEIGHT_GAIN_RESULTS = ClutterResults.modifiedHeights_m;
    const auto MOD_PATH = ClutterResults.modifiedPath;

    const double EXPECTED_HTS = 25.0;
    const double EXPECTED_HRS = 25.0;
    
    EXPECT_NEAR(EXPECTED_HTS,HEIGHT_GAIN_RESULTS.first + MOD_PATH.front().h_asl_m,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_HRS,HEIGHT_GAIN_RESULTS.second+ MOD_PATH.back().h_asl_m,TOLERANCE_STRICT);
}
TEST_F(UrbanProfileTests, calcTimePercentBeta0Test){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EXPECTED_B0 = 7.005407788;
	EXPECT_NEAR(EXPECTED_B0,B0_PERCENT,TOLERANCE_STRICT);
}
TEST_F(UrbanProfileTests, calcFracOverSeaTest){
	// Arrange
    const double VAL_FRAC= K_PATH.calcFracOverSea();
    const double EXPECTED_FRAC = 0.0;
	EXPECT_NEAR(EXPECTED_FRAC,VAL_FRAC,TOLERANCE_STRICT);
}
TEST_F(UrbanProfileTests, calcLongestContiguousInlandDistanceTest){
	// Arrange
    const double VAL_DIST= K_PATH.calcLongestContiguousInlandDistance_km();
    const double EXPECTED_DIST = 5.0;
	EXPECT_NEAR(EXPECTED_DIST,VAL_DIST,TOLERANCE_STRICT);
}

TEST_F(UrbanProfileTests, BasicPropTests_calcPathLossWithGasAndMultipathTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);
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

        EXPECT_NEAR(EXPECTED_LBFSG[freqInd],VAL_LBFSG,TOLERANCE_STRICT);
        EXPECT_NEAR(EXPECTED_LB0P[freqInd],VAL_LB0P,TOLERANCE_STRICT);
        EXPECT_NEAR(EXPECTED_LB0B[freqInd],VAL_LB0B,TOLERANCE_STRICT);
    }
}

TEST_F(UrbanProfileTests, DiffractionLossTests_calcDiffractionLossTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LD50 = {
        0,5.20387644,1.85407598,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    const std::vector<double> EXPECTED_LDP = {
        0,5.20231645,1.85229434,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
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

TEST_F(UrbanProfileTests, TroposcatterLossTests_calcTroposcatterLossTest){
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
TEST_F(UrbanProfileTests, calcHorizonAnglesAndDistances_LineOfSightTest){
    const double FREQ_GHZ = 2.0;
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);

    const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
    const auto mod_path = ClutterResults.modifiedPath;
    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterResults.modifiedHeights_m;

    const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
    const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

    const double EXPECTED_THETA_T = -0.257856304;
    const double EXPECTED_THETA_R = -0.257856304;
    const double EXPECTED_DLT = 2.48;
    const double EXPECTED_DLR = 2.48;

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
TEST_F(UrbanProfileTests, AnomalousProp_calcAnomalousPropLossTest){

    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);
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

