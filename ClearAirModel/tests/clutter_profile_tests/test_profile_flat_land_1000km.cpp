#include "gtest/gtest.h"
#include "test_profile_flat_land_1000km.h"
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
PathProfile::Path FlatLand1000kmProfileTests::K_PATH;

//test height gain model heights
TEST_F(FlatLand1000kmProfileTests, ClutterLossTests_calcHeightGainModelHeightsTest){
    const double FREQ_GHZ = 2;
    const auto ClutterModel = ClearAirModel::ClutterLoss(FREQ_GHZ, K_PATH, 
        HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);

    const auto HEIGHT_GAIN_RESULTS = ClutterModel.getHeightGainModelHeights_m();
    const auto MOD_PATH = ClutterModel.getModifiedPath();
    const double EXPECTED_HTS = 10.0;
    const double EXPECTED_HRS = 10.0;
    EXPECT_NEAR(EXPECTED_HTS,HEIGHT_GAIN_RESULTS.first + MOD_PATH.front().h_asl_m,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_HRS,HEIGHT_GAIN_RESULTS.second+ MOD_PATH.back().h_asl_m,TOLERANCE_STRICT);
}
TEST_F(FlatLand1000kmProfileTests, calcTimePercentBeta0Test){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EXPECTED_B0 = 1.224161171;
	EXPECT_NEAR(EXPECTED_B0,B0_PERCENT,TOLERANCE_STRICT);
}
TEST_F(FlatLand1000kmProfileTests, calcFracOverSeaTest){
	// Arrange
    const double VAL_FRAC= K_PATH.calcFracOverSea();
    const double EXPECTED_FRAC = 0.0;
	EXPECT_NEAR(EXPECTED_FRAC,VAL_FRAC,TOLERANCE_STRICT);
}
TEST_F(FlatLand1000kmProfileTests, calcLongestContiguousInlandDistanceTest){
	// Arrange
    const double VAL_DIST= K_PATH.calcLongestContiguousInlandDistance_km();
    const double EXPECTED_DIST = 1000.0;
	EXPECT_NEAR(EXPECTED_DIST,VAL_DIST,TOLERANCE_STRICT);
}

TEST_F(FlatLand1000kmProfileTests, BasicPropTests_calcPathLossWithGasAndMultipathTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ClearAirModel::ClearAirModelHelpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LBFSG = {
        165.3379426,132.6022376,136.3616371,140.3639758,144.754696,
        149.571306,154.5612894,159.3263438,163.7114343,167.8653361,
        172.0932882,176.8868071,183.3499555,195.5576896,268.9237377,
        277.4859418,358.4305335,574.6456604,165.3379426,165.3379426,
        165.3379426,165.3379426,165.3379426,165.3379426,165.3379426,
        165.3379426,165.3379426,165.3379426,165.3379426,165.3379426,
        165.3379426,165.3379426,165.3379426,165.3379426,165.3379426
    };
    const std::vector<double> EXPECTED_LB0P = {
        165.3165177,132.5808126,136.3402121,140.3425508,144.733271,
        149.549881,154.5398644,159.3049188,163.6900093,167.8439111,
        172.0718632,176.8653821,183.3285305,195.5362646,268.9023128,
        277.4645168,358.4091085,574.6242354,161.1892383,162.6594048,
        163.252878,163.6311321,163.9093697,164.1295714,164.3118188,
        164.4672919,164.6028594,164.7230445,164.8309852,164.9289472,
        165.0186204,165.1012987,165.1779948,165.2495164,165.3165177
    };
    const std::vector<double> EXPECTED_LB0B = {
        161.4037308,128.6680258,132.4274252,136.4297639,140.8204842,
        145.6370942,150.6270776,155.3921319,159.7772225,163.9311242,
        168.1590764,172.9525952,179.4157436,191.6234777,264.9895259,
        273.5517299,354.4963217,570.7114486,161.4037308,161.4037308,
        161.4037308,161.4037308,161.4037308,161.4037308,161.4037308,
        161.4037308,161.4037308,161.4037308,161.4037308,161.4037308,
        161.4037308,161.4037308,161.4037308,161.4037308,161.4037308
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

        //WARNING a higher tolerance is used. 
        EXPECT_NEAR(EXPECTED_LBFSG[freqInd],VAL_LBFSG,TOLERANCE);
        EXPECT_NEAR(EXPECTED_LB0P[freqInd],VAL_LB0P,TOLERANCE);
        EXPECT_NEAR(EXPECTED_LB0B[freqInd],VAL_LB0B,TOLERANCE);
    }
}

TEST_F(FlatLand1000kmProfileTests, DiffractionLossTests_calcDiffractionLossTest){
	// Arrange
    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    const std::vector<double> EXPECTED_LD50 = {
        1048.799584,410.9075607,462.8438751,523.0450456,592.7031384,
        673.1793566,766.0262277,873.0115979,996.1441058,1137.700011,
        1300.252597,1486.709286,1700.369199,1945.177229,2225.804681,
        2547.389506,2915.782884,3047.552268,1048.799584,1048.799584,
        1048.799584,1048.799584,1048.799584,1048.799584,1048.799584,
        1048.799584,1048.799584,1048.799584,1048.799584,1048.799584,
        1048.799584,1048.799584,1048.799584,1048.799584,1048.799584
    };
    const std::vector<double> EXPECTED_LDP = {
        1044.491427,409.362955,461.0656103,520.9994021,590.3515291,
        670.4776369,762.9239131,869.4509839,992.0592778,1133.015735,
        1294.883056,1480.556606,1693.321595,1937.104849,2216.557879,
        2536.798839,2903.654589,3034.874299,661.0582083,747.0243172,
        794.4166973,827.9102982,854.6689385,877.4239106,897.5286242,
        915.7605007,932.6118441,948.4172998,963.4172903,977.7926751,
        991.6851189,1005.209811,1018.463851,1031.53209,1044.491427
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

TEST_F(FlatLand1000kmProfileTests, TroposcatterLossTests_calcTroposcatterLossTest){
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
        272.8392056,229.5307683,235.2393436,241.0359218,247.0642253,
        253.3569498,259.6555114,265.5575686,270.8973845,275.7881778,
        280.4472506,285.1540383,290.4357404,298.0690023,330.8569174,
        340.2669281,410.4203534,611.122357,258.5691372,262.4305231,
        264.1634515,265.3457551,266.2663555,267.0337302,267.701175,
        268.2992565,268.8475015,269.3595329,269.8457003,270.3146643,
        270.7746203,271.2346923,271.707473,272.216377,272.8392056
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
TEST_F(FlatLand1000kmProfileTests, calcHorizonAnglesAndDistances_LineOfSightTest){
    const double FREQ_GHZ = 2.0;
    const double EFF_RADIUS_MED_KM = ClearAirModel::ClearAirModelHelpers::calcMedianEffectiveRadius_km(DN);

    const auto ClutterModel = ClearAirModel::ClutterLoss(FREQ_GHZ, K_PATH, 
            HTG, HRG, TX_CLUTTER_TYPE, RX_CLUTTER_TYPE);
    const auto mod_path = ClutterModel.getModifiedPath();
    const auto [hg_height_tx_m, hg_height_rx_m] = ClutterModel.getHeightGainModelHeights_m();
    const double height_tx_asl_m = hg_height_tx_m + mod_path.front().h_asl_m;
    const double height_rx_asl_m = hg_height_rx_m + mod_path.back().h_asl_m;

    const double EXPECTED_THETA_T = -1.442104943;
    const double EXPECTED_THETA_R = -1.442104943;
    const double EXPECTED_DLT = 14.0;
    const double EXPECTED_DLR = 14.0;

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
TEST_F(FlatLand1000kmProfileTests, AnomolousProp_calcAnomolousPropLossTest){

    const double B0_PERCENT = K_PATH.calcTimePercentBeta0(INPUT_LAT);
    const double EFF_RADIUS_MED_KM = ClearAirModel::ClearAirModelHelpers::calcMedianEffectiveRadius_km(DN);
    const double SEA_FRAC = K_PATH.calcFracOverSea();

    //basic transmission loss from ducting/layer refraction
    const std::vector<double> EXPECTED_LBA = {
        334.9166351,296.0965295,297.4275348,297.4943355,296.6048981,
        296.6457139,307.2425464,318.4258409,330.1579823,342.7221578,
        356.5774701,372.3915652,391.4701242,418.1188986,508.0158111,
        535.5011311,638.1072832,862.0676781,307.5854244,316.7485585,
        320.6190061,323.1400854,325.0222282,326.5287336,327.7870663,
        328.868861,329.8184702,330.6652987,331.4298623,332.1270536,
        332.7680309,333.3613749,333.9138279,334.4307857,334.9166351
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

TEST_F(FlatLand1000kmProfileTests, calcP452TotalAttenuationTest){
    //basic loss prediction not exceeded for p percent
    const std::vector<double> EXPECTED_LOSS = {
        272.8392056,229.5307683,235.2393436,241.0359218,247.0642253,
        253.3569498,259.6555114,265.5575686,270.8973845,275.7881778,
        280.4472506,285.1540383,290.4357404,298.0690023,330.8569174,
        340.2669281,410.4203534,611.122357,258.5691372,262.4305231,
        264.1634515,265.3457551,266.2663555,267.0337302,267.701175,
        268.2992565,268.8475015,269.3595329,269.8457003,270.3146643,
        270.7746203,271.2346923,271.707473,272.216377,272.8392056
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

