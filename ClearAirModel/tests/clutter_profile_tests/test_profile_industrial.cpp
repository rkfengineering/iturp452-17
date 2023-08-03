#include "gtest/gtest.h"
#include "test_profile_industrial.h"
#include "ClearAirModel/ClearAirModelHelpers.h"
#include "ClearAirModel/BasicProp.h"
#include "ClearAirModel/DiffractionLoss.h"
#include "ClearAirModel/TropoScatter.h"
#include "ClearAirModel/AnomalousProp.h"
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
PathProfile::Path IndustrialProfileTests::K_PATH;

//test height gain model heights
TEST_F(IndustrialProfileTests, ClutterLossTests_calcHeightGainModelHeightsTest){
    const double FREQ_GHZ = 2;
    const auto ClutterModel = ClearAirModel::ClutterLoss(FREQ_GHZ, K_PATH, 
        HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);

    const auto HEIGHT_GAIN_RESULTS = ClutterModel.getHeightGainModelHeights_m();
    const auto MOD_PATH = ClutterModel.getModifiedPath();
    const double EXPECTED_HTS = 20.0;
    const double EXPECTED_HRS = 20.0;
    EXPECT_NEAR(EXPECTED_HTS,HEIGHT_GAIN_RESULTS.first + MOD_PATH.front().h_asl_m,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_HRS,HEIGHT_GAIN_RESULTS.second+ MOD_PATH.back().h_asl_m,TOLERANCE_STRICT);
}
TEST_F(IndustrialProfileTests, calcTimePercentBeta0Test){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EXPECTED_B0 = 7.005407788;
	EXPECT_NEAR(EXPECTED_B0,B0_PERCENT,TOLERANCE_STRICT);
}
TEST_F(IndustrialProfileTests, calcFracOverSeaTest){
	// Arrange
    const double VAL_FRAC= K_PATH.calcFracOverSea();
    const double EXPECTED_FRAC = 0.0;
	EXPECT_NEAR(EXPECTED_FRAC,VAL_FRAC,TOLERANCE_STRICT);
}
TEST_F(IndustrialProfileTests, calcLongestContiguousInlandDistanceTest){
	// Arrange
    const double VAL_DIST= K_PATH.calcLongestContiguousInlandDistance_km();
    const double EXPECTED_DIST = 5.0;
	EXPECT_NEAR(EXPECTED_DIST,VAL_DIST,TOLERANCE_STRICT);
}

TEST_F(IndustrialProfileTests, BasicPropTests_calcPathLossWithGasAndMultipathTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ClearAirModel::ClearAirModelHelpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LBFSG = {
        112.2584165,86.20491256,89.72790186,93.25208156,96.77816432,
        100.306334,103.8353531,107.3632701,110.8893253,114.4142477,
        117.9395329,121.4675893,125.003827,128.5682131,132.432275,
        135.978798,139.8799948,142.0858264,112.2584165,112.2584165,
        112.2584165,112.2584165,112.2584165,112.2584165,112.2584165,
        112.2584165,112.2584165,112.2584165,112.2584165,112.2584165,
        112.2584165,112.2584165,112.2584165,112.2584165,112.2584165
    };
    const std::vector<double> EXPECTED_LB0P = {
        112.2495797,86.19607572,89.71906501,93.24324471,96.76932748,
        100.2974971,103.8265163,107.3544333,110.8804885,114.4054108,
        117.930696,121.4587525,124.9949902,128.5593763,132.4234382,
        135.9699612,139.8711579,142.0769895,110.5472625,111.1536401,
        111.3984211,111.5544339,111.6691944,111.7600177,111.8351866,
        111.8993122,111.9552277,112.0047987,112.0493194,112.0897243,
        112.1267104,112.1608115,112.1924452,112.2219446,112.2495797
    };
    const std::vector<double> EXPECTED_LB0B = {
        111.3987589,85.34525495,88.86824425,92.39242394,95.91850671,
        99.44667634,102.9756955,106.5036125,110.0296677,113.55459,
        117.0798752,120.6079317,124.1441694,127.7085555,131.5726174,
        135.1191404,139.0203372,141.2261688,111.3987589,111.3987589,
        111.3987589,111.3987589,111.3987589,111.3987589,111.3987589,
        111.3987589,111.3987589,111.3987589,111.3987589,111.3987589,
        111.3987589,111.3987589,111.3987589,111.3987589,111.3987589
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

TEST_F(IndustrialProfileTests, DiffractionLossTests_calcDiffractionLossTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LD50 = {
        0,8.68406638,5.45980554,2.13819017,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    const std::vector<double> EXPECTED_LDP = {
        0,8.68240951,5.4579253,2.13604177,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
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

TEST_F(IndustrialProfileTests, TroposcatterLossTests_calcTroposcatterLossTest){
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
        161.9975672,125.2080198,130.6794488,135.997028,141.1614617,
        146.17291,151.0301061,155.7310783,160.2750146,164.6624697,
        168.8945084,172.9724999,176.8990276,180.6827969,184.4355431,
        187.9194567,191.5467319,193.5707559,147.7274988,151.5888847,
        153.321813,154.5041167,155.4247171,156.1920918,156.8595366,
        157.4576181,158.0058631,158.5178945,159.0040619,159.4730259,
        159.9329819,160.3930539,160.8658346,161.3747386,161.9975672
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
TEST_F(IndustrialProfileTests, calcHorizonAnglesAndDistances_LineOfSightTest){
    const double FREQ_GHZ = 2.0;
    const double EFF_RADIUS_MED_KM = ClearAirModel::ClearAirModelHelpers::calcMedianEffectiveRadius_km(DN);

    const auto ClutterModel = ClearAirModel::ClutterLoss(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
    const auto mod_path = ClutterModel.getModifiedPath();
    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterModel.getHeightGainModelHeights_m();
    const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
    const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

    const double EXPECTED_THETA_T = -0.254737074;
    const double EXPECTED_THETA_R = -0.254737074;
    const double EXPECTED_DLT = 2.45;
    const double EXPECTED_DLR = 2.45;

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
TEST_F(IndustrialProfileTests, AnomalousProp_calcAnomalousPropLossTest){

    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ClearAirModel::ClearAirModelHelpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    //basic transmission loss from ducting/layer refraction
    const std::vector<double> EXPECTED_LBA = {
        183.9648729,190.511369,188.3406083,184.1913505,178.1584489,
        172.0127904,175.5418095,179.0697265,182.5957817,186.1207041,
        189.6459893,193.1740457,196.7102834,200.2746695,204.1387314,
        207.6852545,211.5864512,213.7922828,111.2584583,117.1554841,
        122.2994335,127.1787468,131.9029892,136.5195044,141.0542793,
        145.523581,149.9384838,154.3069816,158.6350985,162.9275271,
        167.1880201,171.4196431,175.6249446,179.8060753,183.9648729
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
        const auto AnomalousModel = ClearAirModel::AnomalousProp(
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
        const double LOSS_VAL = AnomalousModel.getAnomalousPropLoss_dB();
        EXPECT_NEAR(EXPECTED_LBA[freqInd],LOSS_VAL,TOLERANCE);
    }
}

TEST_F(IndustrialProfileTests, calcP452TotalAttenuationTest){

    //basic loss prediction not exceeded for p percent
    const std::vector<double> EXPECTED_LOSS = {
        143.4695904,102.2449514,102.6093824,103.0674969,105.9967285,
        120.1204379,134.5676804,138.5728015,142.1004989,145.6254216,
        149.1507068,152.6787633,156.2150009,159.7793871,163.643449,
        167.189972,171.0911687,173.2970003,141.7729993,142.3793769,
        142.6241579,142.7374797,142.8279321,142.9048509,142.972811,
        143.0344404,143.0914031,143.1448304,143.195535,143.2441282,
        143.2910889,143.3368065,143.3816092,143.4257839,143.4695904
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

