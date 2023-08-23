#include "gtest/gtest.h"
#include "land_70km_testclass.h"
#include "MainModel/Helpers.h"
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
PathProfile::Path Land70kmProfileTests::K_PATH;

//test height gain model heights
TEST_F(Land70kmProfileTests, ClutterModelTests_calcHeightGainModelHeightsTest){
    const double FREQ_GHZ = 2;
    const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
    const auto HEIGHT_GAIN_RESULTS = ClutterResults.modifiedHeights_m;
    const auto MOD_PATH = ClutterResults.modifiedPath;
    const double EXPECTED_HTS = 837.0;
    const double EXPECTED_HRS = 702.0;
    EXPECT_NEAR(EXPECTED_HTS,HEIGHT_GAIN_RESULTS.first + MOD_PATH.front().h_asl_m,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_HRS,HEIGHT_GAIN_RESULTS.second+ MOD_PATH.back().h_asl_m,TOLERANCE_STRICT);
}
TEST_F(Land70kmProfileTests, calcTimePercentBeta0Test){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EXPECTED_B0 = 2.563929559;
	EXPECT_NEAR(EXPECTED_B0,B0_PERCENT,TOLERANCE_STRICT);
}
TEST_F(Land70kmProfileTests, calcFracOverSeaTest){
	// Arrange
    const double VAL_FRAC= K_PATH.calcFracOverSea();
    const double EXPECTED_FRAC = 0.0;
	EXPECT_NEAR(EXPECTED_FRAC,VAL_FRAC,TOLERANCE_STRICT);
}
TEST_F(Land70kmProfileTests, calcLongestContiguousInlandDistanceTest){
	// Arrange
    const double VAL_DIST= K_PATH.calcLongestContiguousInlandDistance_km();
    const double EXPECTED_DIST = 69.94042916;
	EXPECT_NEAR(EXPECTED_DIST,VAL_DIST,TOLERANCE_STRICT);
}

TEST_F(Land70kmProfileTests, BasicPropTests_calcPathLossWithGasAndMultipathTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LBFSG = {
        135.7989848,109.3087267,112.8471679,116.4026005,119.9851967,
        123.5975798,127.2220888,130.8308661,134.4130686,137.9791016,
        141.5503136,145.1610818,148.8886248,153.0179473,161.4247165,
        165.2990686,174.2358854,190.4295416,135.7989848,135.7989848,
        135.7989848,135.7989848,135.7989848,135.7989848,135.7989848,
        135.7989848,135.7989848,135.7989848,135.7989848,135.7989848,
        135.7989848,135.7989848,135.7989848,135.7989848,135.7989848
    };
    const std::vector<double> EXPECTED_LB0P = {
        134.6229822,108.1327241,111.6711653,115.2265979,118.8091941,
        122.4215772,126.0460862,129.6548636,133.237066,136.803099,
        140.3743111,143.9850792,147.7126222,151.8419447,160.2487139,
        164.123066,173.0598828,189.253539,132.9405029,133.9534564,
        134.3623629,134.6229822,134.8146895,134.9664099,135.0919794,
        135.1991013,135.292508,135.3753163,135.4496881,135.5171845,
        135.5789698,135.6359357,135.6887798,135.7380585,135.784222
    };
    const std::vector<double> EXPECTED_LB0B = {
        133.6284789,107.1382208,110.6766621,114.2320946,117.8146908,
        121.427074,125.051583,128.6603603,132.2425627,135.8085957,
        139.3798078,142.9905759,146.7181189,150.8474414,159.2542106,
        163.1285627,172.0653795,188.2590357,133.6284789,133.6284789,
        133.6284789,133.6284789,133.6284789,133.6284789,133.6284789,
        133.6284789,133.6284789,133.6284789,133.6284789,133.6284789,
        133.6284789,133.6284789,133.6284789,133.6284789,133.6284789
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

TEST_F(Land70kmProfileTests, DiffractionLossTests_calcDiffractionLossTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LD50 = {
        58.42626086,50.9239042,50.67353983,50.8186498,51.37763005,
        52.35766638,53.74409809,55.49162017,57.52877055,59.97800622,
        63.08020767,66.91393844,71.53688544,77.02634313,83.4765291,
        90.99873253,99.72296601,102.8624531,58.42626086,58.42626086,
        58.42626086,58.42626086,58.42626086,58.42626086,58.42626086,
        58.42626086,58.42626086,58.42626086,58.42626086,58.42626086,
        58.42626086,58.42626086,58.42626086,58.42626086,58.42626086
    };
    const std::vector<double> EXPECTED_LDP = {
        51.13018607,48.43167876,47.90764595,47.65904647,47.69842342,
        48.01287719,48.63400439,49.54567949,50.65560047,51.92352028,
        53.46073197,55.40331252,57.74302964,60.47643212,63.61354332,
        67.17509053,71.20974707,72.66508627,47.3285835,48.45848488,
        50.02387701,51.13018607,52.01403627,52.76564347,53.42971142,
        54.03191869,54.58852632,55.1105878,55.60604439,56.08086996,
        56.53974378,56.98647057,57.42425764,57.85590758,58.28396044
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

TEST_F(Land70kmProfileTests, TroposcatterLossTests_calcTroposcatterLossTest){
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
        196.4720423,159.2564157,164.7433447,170.0922316,175.3131284,
        180.4083219,185.3597266,190.1392026,194.7351581,199.1555117,
        203.4154617,207.5345518,211.5496546,215.5850295,221.2355345,
        225.1067819,233.0822527,248.0920078,189.6954243,193.5568103,
        195.2897386,196.4720423,197.3926427,198.1600174,198.8274622,
        199.4255437,199.9737887,200.4858201,200.9719875,201.4409515,
        201.9009075,202.3609795,202.8337602,203.3426642,203.9654928
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
TEST_F(Land70kmProfileTests, calcHorizonAnglesAndDistances_LineOfSightTest){
    const double FREQ_GHZ = 2.0;
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);

    const auto ClutterResults = ClutterModel::calculateClutterModel(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
        const auto mod_path = ClutterResults.modifiedPath;
    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterResults.modifiedHeights_m;

    const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
    const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

    const double EXPECTED_THETA_T = 0.698535227;
    const double EXPECTED_THETA_R = 16.76431411;
    const double EXPECTED_DLT = 9.227522888;
    const double EXPECTED_DLR = 1.188393099;

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
TEST_F(Land70kmProfileTests, AnomalousProp_calcAnomalousPropLossTest){

    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ITUR_P452::Helpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    //basic transmission loss from ducting/layer refraction
    const std::vector<double> EXPECTED_LBA = {
        195.0165905,184.0632766,183.6895346,181.50586,177.6215949,
        173.8164872,179.8699128,186.0910181,192.4827518,199.072745,
        205.9034773,213.0348512,220.5745335,228.8434869,241.7595081,
        250.5618318,264.9022351,282.9774689,163.0857196,178.1857615,
        187.596647,195.0165905,201.3274366,206.9071539,211.9587675,
        216.6061296,220.9313249,224.9920426,228.8306429,232.4793066,
        235.9631559,239.3022444,242.5128801,245.6085351,248.6004903
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

TEST_F(Land70kmProfileTests, calcP452TotalAttenuationTest){

    //basic loss prediction not exceeded for p percent
    const std::vector<double> EXPECTED_LOSS = {
        185.7376292,156.0123535,159.3863008,162.808431,166.4702998,
        170.4125873,174.6642693,179.186495,183.8779844,188.7088692,
        193.8088582,199.3379869,205.328257,211.8827332,220.6686567,
        224.9848195,233.0697205,248.0882834,163.0953075,178.186287,
        184.3719632,185.7376292,186.8120418,187.7142957,188.5029097,
        189.2112561,189.8603267,190.4642964,191.0332804,191.5748351,
        192.0948396,192.5980505,193.0884817,193.569723,194.0455355
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        const auto p452Model = ITUR_P452::TotalClearAirAttenuation(
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
        double LOSS_VAL = p452Model.calcTotalClearAirAttenuation();
        EXPECT_NEAR(EXPECTED_LOSS[freqInd],LOSS_VAL,TOLERANCE);
    }
}

