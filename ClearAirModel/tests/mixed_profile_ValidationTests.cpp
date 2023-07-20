#include "gtest/gtest.h"
#include "ClearAirModel/EffectiveEarth.h"
#include "ClearAirModel/BasicProp.h"
#include "ClearAirModel/DiffractionLoss.h"
#include <filesystem>

//Validation data from ITU validation spreadsheets in results folder
//in ITU validation dataset titled "Validation examples for software implementation of Recommendation ITU-R P.452"
//https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/R19-WP3M-C-0364!N18-P2!ZIP-E.zip

//Intermediate results from the transmission losses tab of the spreadsheets in the validation examples folder are inaccurate.
//They do not account for the 0.1dB difference in the free space path loss constant changed between 452-16 and 452-17. 
//As such, results from the spreadsheets in the results folder are preferred


//If there's a way to use relative paths that would be nice. The paths from the validation data have been isolated as separate
//csv files and put in a subdirectory in the tests folder
const std::filesystem::path testPathFileDir("/home/ayeh/itu/ituModels/iturp452/ClearAirModel/tests/test_paths");

namespace {
	// Use when expected an exact match
	double constexpr TOLERANCE_STRICT = 1.0e-6;
    // Use when different speed of light constant used
    // The validation data will match to a stricter tolerance if 2.998*1e8 is used for speed of light
    double constexpr TOLERANCE = 1.0e-3;
}

//From intermediate values in test_result_mixed_109.csv
TEST(EffectiveEarthTests, EffectiveEarthTests_effRadiusp50km){
    const double INPUT_DN = 53;
    const double EXPECTED_AE = 9617.759615;
	const double VAL_AE = EffectiveEarth::calcMedianEffectiveRadius_km(INPUT_DN);

	EXPECT_NEAR(EXPECTED_AE,VAL_AE,TOLERANCE_STRICT);
}

//check path fraction calculation (omega) using mixed terrain path
TEST(ProfilePathTests, ProfilePath_calcFracOverSea){
    const PathProfile::Path p(testPathFileDir/std::filesystem::path("test_profile_mixed_109km.csv"));
    const double VAL_FRAC= p.calcFracOverSea();
    const double EXPECTED_FRAC = 0.394495413;
	EXPECT_NEAR(EXPECTED_FRAC,VAL_FRAC,TOLERANCE_STRICT);
}

//check beta0 using mixed terrain path
TEST(ProfilePathTests, ProfilePath_calcTimePercentBeta0){
    const PathProfile::Path p(testPathFileDir/std::filesystem::path("test_profile_mixed_109km.csv"));
    const double PHI_T = 51.2;
    const double PHI_R = 50.73;
    const double INPUT_LAT = (PHI_T+PHI_R)/2.0;

    const double VAL_B0= p.calcTimePercentBeta0(INPUT_LAT);
    const double EXPECTED_B0 = 3.282731389;
	EXPECT_NEAR(EXPECTED_B0,VAL_B0,TOLERANCE_STRICT);
}

//check horizon elevation angle and distance calculations using mixed terrain path
//transhorizon path tested
TEST(EffectiveEarthTests, EffectiveEarthTests_HorizonAnglesAndDistances_Transhorizon){
    const PathProfile::Path p(testPathFileDir/std::filesystem::path("test_profile_mixed_109km.csv"));
    const double HTG = 10;
    const double HRG = 10;
    const double DN = 53;
    const double FREQ_GHZ = 0.2;

    const double HTS_MASL = HTG+p.front().h_asl_m;
    const double HRS_MASL = HRG+p.back().h_asl_m;
    const double EFF_RADIUS_MED_KM = EffectiveEarth::calcMedianEffectiveRadius_km(DN);

    const double EXPECTED_THETA_T = -0.6342118;
    const double EXPECTED_THETA_R = -1.390039674;
    const double EXPECTED_DLT = 28;
    const double EXPECTED_DLR = 11;

    const auto HORIZON_VALS = EffectiveEarth::calcHorizonAnglesAndDistances(
        p,
        HTS_MASL,
        HRS_MASL,
        EFF_RADIUS_MED_KM,
        FREQ_GHZ
    );
    EXPECT_NEAR(EXPECTED_THETA_T,HORIZON_VALS.horizonElevation_mrad.tx_val,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_THETA_R,HORIZON_VALS.horizonElevation_mrad.rx_val,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_DLT,HORIZON_VALS.horizonDist_km.tx_val,TOLERANCE_STRICT);
    EXPECT_NEAR(EXPECTED_DLR,HORIZON_VALS.horizonDist_km.rx_val,TOLERANCE_STRICT);
}

TEST(BasicPropTests, BasicPropTests_calcPathLossWithGasAndMultipath_dB){
	// Arrange
	const std::vector<double> FREQ_GHZ_LIST = {
		0.2,0.1,0.15,0.225,0.3375,0.50625,0.759375,1.1390625,
		1.70859375,2.562890625,3.844335938,5.766503906,8.649755859,
		12.97463379,19.46195068,29.19292603,43.78938904,50,0.2,0.2,
		0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2
	};
    const std::vector<double> P_LIST = {
        0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
        0.1,0.1,0.1,1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,
    };

    const PathProfile::Path p(testPathFileDir/std::filesystem::path("test_profile_mixed_109km.csv"));
    const double HTG = 10;
    const double HRG = 10;
    const double TEMP_C = 15;
    const double DRY_PRESSURE_HPA = 1013;
    const double DN = 53;
    const double PHI_T = 51.2;
    const double PHI_R = 50.73;

    const double HTS_MASL = HTG+p.front().h_asl_m;
    const double HRS_MASL = HRG+p.back().h_asl_m;
    const double TEMP_K = TEMP_C + 273.15;
    const double SEA_FRAC = p.calcFracOverSea();
    const double EFF_RADIUS_MED_KM = EffectiveEarth::calcMedianEffectiveRadius_km(DN);
    const double INPUT_LAT = (PHI_T+PHI_R)/2.0;
    const double VAL_B0= p.calcTimePercentBeta0(INPUT_LAT);

    const std::vector<double> EXPECTED_LB0P = {
        112.3752248,106.2952854,109.8429897,113.4171732,117.0337313,
        120.6968708,124.3792214,128.0374844,131.6550584,135.2491031,
        138.8552054,142.5325022,146.4162993,151.0088874,163.1280578,
        167.2307246,179.519313,204.526306,114.9225958,116.456266,117.0753741,
        117.4699669,117.7602229,117.9899362,118.1800557,118.3422445,118.4836677,
        118.6090443,118.7216474,118.8238409,118.9173875,119.0036371,119.0836461,
        119.158257,119.2281524
    };
    const std::vector<double> EXPECTED_LB0B = {
        116.2376388,110.1576994,113.7054037,117.2795872,120.8961453,124.5592848,
        128.2416354,131.8998984,135.5174724,139.1115171,142.7176194,146.3949162,
        150.2787133,154.8713014,166.9904718,171.0931386,183.381727,208.38872,
        116.2376388,116.2376388,116.2376388,116.2376388,116.2376388,116.2376388,
        116.2376388,116.2376388,116.2376388,116.2376388,116.2376388,116.2376388,
        116.2376388,116.2376388,116.2376388,116.2376388,116.2376388
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {

        const auto HORIZON_VALS = EffectiveEarth::calcHorizonAnglesAndDistances(
            p, HTS_MASL, HRS_MASL, EFF_RADIUS_MED_KM, FREQ_GHZ_LIST[freqInd]
        );
        const double VAL_LB0P = BasicProp::calcPathLossWithGasAndMultipath_dB(
            p.back().d_km,
            HTS_MASL,
            HRS_MASL,
            FREQ_GHZ_LIST[freqInd],
            TEMP_K,
            DRY_PRESSURE_HPA,
            SEA_FRAC,
            HORIZON_VALS.horizonDist_km.tx_val,
            HORIZON_VALS.horizonDist_km.rx_val,
            P_LIST[freqInd]
        );
        const double VAL_LB0B = BasicProp::calcPathLossWithGasAndMultipath_dB(
            p.back().d_km,
            HTS_MASL,
            HRS_MASL,
            FREQ_GHZ_LIST[freqInd],
            TEMP_K,
            DRY_PRESSURE_HPA,
            SEA_FRAC,
            HORIZON_VALS.horizonDist_km.tx_val,
            HORIZON_VALS.horizonDist_km.rx_val,
            VAL_B0
        );
        EXPECT_NEAR(EXPECTED_LB0P[freqInd],VAL_LB0P,TOLERANCE_STRICT);
        EXPECT_NEAR(EXPECTED_LB0B[freqInd],VAL_LB0B,TOLERANCE_STRICT);
    }
}

TEST(DiffractionLossTests, DiffractionLoss_sphericalEarthDiffLossTest){
	// Arrange
	const std::vector<double> FREQ_GHZ_LIST = {
		0.2,0.1,0.15,0.225,0.3375,0.50625,0.759375,1.1390625,
		1.70859375,2.562890625,3.844335938,5.766503906,8.649755859,
		12.97463379,19.46195068,29.19292603,43.78938904,50,0.2,0.2,
		0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2
	};

    const PathProfile::Path p(testPathFileDir/std::filesystem::path("test_profile_mixed_109km.csv"));
    const double HTG = 10;
    const double HRG = 10;
    const double DN = 53;
    const auto POL = Enumerations::PolarizationType::HorizontalPolarized;

    const double HTS_MASL = HTG+p.front().h_asl_m;
    const double HRS_MASL = HRG+p.back().h_asl_m;
    const double SEA_FRAC = p.calcFracOverSea();
    const double EFF_RADIUS_MED_KM = EffectiveEarth::calcMedianEffectiveRadius_km(DN);

    const std::vector<double> EXPECTED_LDSPH = {
        32.52266992,32.18219728,32.23540955,32.71355093,33.68370872,
        35.22094531,37.28531414,39.86569196,43.08667459,47.04755022,
        51.79540079,57.39874993,63.94574231,71.54393877,80.32127304,
        90.4278295,102.0382751,106.1957983,32.52266992,32.52266992,
        32.52266992,32.52266992,32.52266992,32.52266992,32.52266992,
        32.52266992,32.52266992,32.52266992,32.52266992,32.52266992,
        32.52266992,32.52266992,32.52266992,32.52266992,32.52266992
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        //input effective antanna heights relative to ground, frequency in GHz
        const EffectiveEarth::TxRxPair height_eff_m = EffectiveEarth::calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m(p,HTG,HRG);
        
        const double LDSPH = DiffractionLoss::calcSphericalEarthDiffractionLoss_dB(
            p.back().d_km, 
            HTS_MASL-height_eff_m.tx_val, //effective height relative to ground
            HRS_MASL-height_eff_m.rx_val,
            EFF_RADIUS_MED_KM, 
            FREQ_GHZ_LIST[freqInd],
            SEA_FRAC,
            POL
        );

        EXPECT_NEAR(EXPECTED_LDSPH[freqInd],LDSPH,TOLERANCE_STRICT);
    }
}

TEST(DiffractionLossTests, DiffractionLoss_calcDiffractionLoss_dBTest){
	// Arrange
	const std::vector<double> FREQ_GHZ_LIST = {
		0.2,0.1,0.15,0.225,0.3375,0.50625,0.759375,1.1390625,
		1.70859375,2.562890625,3.844335938,5.766503906,8.649755859,
		12.97463379,19.46195068,29.19292603,43.78938904,50,0.2,0.2,
		0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2
	};
    const std::vector<double> P_LIST = {
        0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
        0.1,0.1,0.1,1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,
    };

    const PathProfile::Path p(testPathFileDir/std::filesystem::path("test_profile_mixed_109km.csv"));
    const double HTG = 10;
    const double HRG = 10;
    const double DN = 53;
    const double PHI_T = 51.2;
    const double PHI_R = 50.73;
    const auto POL = Enumerations::PolarizationType::HorizontalPolarized;

    const double HTS_MASL = HTG+p.front().h_asl_m;
    const double HRS_MASL = HRG+p.back().h_asl_m;
    const double SEA_FRAC = p.calcFracOverSea();
    const double INPUT_LAT = (PHI_T+PHI_R)/2.0;
    const double VAL_B0= p.calcTimePercentBeta0(INPUT_LAT);

    const std::vector<double> EXPECTED_LD50 = {
        40.80389052,39.42467495,40.11782021,41.14276783,42.54968562,44.41049424,
        46.69475124,49.40947074,52.70077117,56.68831231,61.4357559,67.02368092,
        73.54824404,81.12192602,89.87544161,99.96027974,111.5516773,115.7035738,
        40.80389052,40.80389052,40.80389052,40.80389052,40.80389052,40.80389052,
        40.80389052,40.80389052,40.80389052,40.80389052,40.80389052,40.80389052,
        40.80389052,40.80389052,40.80389052,40.80389052,40.80389052
    };
    const std::vector<double> EXPECTED_LDP = {
        29.85687048,30.54877594,30.01754059,29.8352037,29.9407675,30.36818397,
        31.08570466,31.98494816,33.00678387,34.44136384,36.24781304,38.0529047,
        39.85412486,41.65042978,43.44163699,45.22802743,47.01010003,47.59219178,
        29.85687048,30.39262859,32.0276681,33.18319902,34.10637345,34.89142114,
        35.58503477,36.21403543,36.79540763,37.34069664,37.85819707,38.35414857,
        38.83343858,39.30004113,39.7573062,40.2081611,40.65525888
    };

    for (uint32_t freqInd = 0; freqInd < FREQ_GHZ_LIST.size(); freqInd++) {
        //input effective antanna heights relative to ground, frequency in GHz
        const EffectiveEarth::TxRxPair height_eff_m = EffectiveEarth::calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m(p,HTG,HRG);
        
        const auto LOSS_VAL = DiffractionLoss::calcDiffractionLoss_dB(
            p, 
            HTS_MASL,
            HRS_MASL,
            height_eff_m.tx_val,
            height_eff_m.rx_val,
            FREQ_GHZ_LIST[freqInd],
            SEA_FRAC,
            P_LIST[freqInd],
            VAL_B0,
            DN,
            POL
        );

        EXPECT_NEAR(EXPECTED_LD50[freqInd],LOSS_VAL.diff_loss_p50_dB,TOLERANCE);
        EXPECT_NEAR(EXPECTED_LDP[freqInd],LOSS_VAL.diff_loss_p_dB,TOLERANCE);
    }
}