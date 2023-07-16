#include "gtest/gtest.h"
#include "DiffractionLoss.h"
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
    //assumption for validation data
    double constexpr ae = 8500; //km
}

//TODO move this stuff to a test case set up class
static const std::vector<PathProfile::Path> PROFILE_LIST = {
    PathProfile::Path(testPathFileDir/std::filesystem::path("path1.csv")),
    PathProfile::Path(testPathFileDir/std::filesystem::path("path2.csv")),
    PathProfile::Path(testPathFileDir/std::filesystem::path("path3.csv")),
    PathProfile::Path(testPathFileDir/std::filesystem::path("path4.csv")),
};

TEST(EffectiveEarthTests, EffectiveEarthTests_diffractionModelHeights){
    // Arrange
	const std::vector<int> PATH_LIST = {
		1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4
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
        const EffectiveEarth::HeightPair EFF_HEIGHT = EffectiveEarth::smoothEarthHeights_diffractionModel(
            p,HTS_LIST[pathInd],HRS_LIST[pathInd]);
        const double HTS_AMSL = HTS_LIST[pathInd]+p.front().h_masl;
        const double HRS_AMSL = HRS_LIST[pathInd]+p.back().h_masl;
        //Equation 38 to calculate modified antenna heights
        EXPECT_NEAR(EXPECTED_HTSP[pathInd],HTS_AMSL-EFF_HEIGHT.tx_val,TOLERANCE);
        EXPECT_NEAR(EXPECTED_HRSP[pathInd],HRS_AMSL-EFF_HEIGHT.rx_val,TOLERANCE);
    }
}

TEST(DiffractionLossTests, DiffractionLoss_bullingtonTest){
    // Arrange
	const std::vector<int> PATH_LIST = {
		1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4
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
    
    const std::vector<double> EXPECTED_LBA = {
        32.84284766, 37.229036, 31.22595813, 24.87792613, 24.78223928,
        36.20419305,44.81240928,34.43917939,28.37368149,34.61746032,
        17.62040631,25.1016907,20.25488733,10.1899085,16.27593106,
        7.8709851352,0,15.30876397,6.526730443,0
    };

    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); ++pathInd) {
        //hts,hrs needs to be in m asl, frequency in GHz
        const PathProfile::Path p = PROFILE_LIST[PATH_LIST[pathInd]-1];
        const double LBA = DiffractionLoss::bullLoss(p,
            HTS_LIST[pathInd]+p.front().h_masl,
            HRS_LIST[pathInd]+p.back().h_masl,
            ae, FREQ_MHZ_LIST[pathInd]/1000.0);

        //suggested tolerance 0.1 dB from Notes page of spreadsheet. 
        //Validation data uses 3e8 for speed of light, model uses 2.998e8
        EXPECT_NEAR(EXPECTED_LBA[pathInd],LBA,0.1);
    }
}

//Not sure what polarization type to use
//not specified by validation data
//Assume fraction over sea = 0 since the profile heights in amsl are positive values

//Two tests are excluded (row 6 and row 17) since they take a different branch in the code
//for the spherical earth calculations, resulting in a different modified effective earth radius being used
//for the intermediate calculations that produce the values used to validate this first term function
TEST(DiffractionLossTests, DiffractionLoss_sphericalEarthFTTest){
    // Arrange
	const std::vector<int> PATH_LIST = {
		1,1,1,1,2,2,2,2,2,
	};
    const std::vector<double> HTS_LIST = {
        30,50,20,70,
        30,50,20,40,70
    };
    const std::vector<double> HRS_LIST = {
        30,10,20,5,
        30,10,20,50,5
    };
    const std::vector<double> FREQ_MHZ_LIST = {
        1000,2500,600,150,
        1000,2500,600,200,150
    };

    const std::vector<double> EXPECTED_FX = {
        -59.40894528,-85.5189285,-48.12283666,-26.15770548,
        -102.9062648,-145.2758435,-84.49345195,-54.44670331,-48.32185499
    };
    const std::vector<double> EXPECTED_GY1 = {
        39.52373705,62.46549974,30.34301718,16.59625588,
        7.40810458,23.83868156,0.127517272,-1.774893492,1.161982247
    };
    const std::vector<double> EXPECTED_GY2 = {
        4.536624536,-0.643878137,-3.134029165,-23.57738733,
        26.04381871,35.89852443,18.04945817,10.466862,3.806254354
    };
    
    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
        //input effective antanna heights relative to ground, frequency in GHz
        const PathProfile::Path p = PROFILE_LIST[PATH_LIST[pathInd]-1];
        const EffectiveEarth::HeightPair height_eff_m = EffectiveEarth::smoothEarthHeights_diffractionModel(
            p,HTS_LIST[pathInd],HRS_LIST[pathInd]);
        
        const double FT = DiffractionLoss::se_first_term(p.back().d_km-p.front().d_km, 
            HTS_LIST[pathInd]+p.front().h_masl-height_eff_m.tx_val, //effective height relative to ground
            HRS_LIST[pathInd]+p.back().h_masl-height_eff_m.rx_val,
            ae, FREQ_MHZ_LIST[pathInd]/1000.0,
            0,Enumerations::PolarizationType::HorizontalPolarized); //Assume land, horizontal pol
        const double EXPECTED_FT = -EXPECTED_FX[pathInd] - EXPECTED_GY1[pathInd] - EXPECTED_GY2[pathInd];

        //suggested tolerance 0.1 dB from Notes page of spreadsheet. 
        //Validation data uses 3e8 for speed of light, model uses 2.998e8
        EXPECT_NEAR(EXPECTED_FT,FT,0.1);
    }
}

TEST(DiffractionLossTests, DiffractionLoss_sphericalEarthTest){
    // Arrange
	const std::vector<int> PATH_LIST = {
		1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4
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

    const std::vector<double> EXPECTED_LSPH = {
        15.3485837,23.6973069,20.91384864,14.32933641,33.13883693,
        69.45434151,85.53863749,66.31647651,45.7547348,43.35361839,
        0,0,0,0,10.63765222,
        0,0,0,0,0
    };

    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
        //input effective antanna heights relative to ground, frequency in GHz
        const PathProfile::Path p = PROFILE_LIST[PATH_LIST[pathInd]-1];
        const EffectiveEarth::HeightPair height_eff_m = EffectiveEarth::smoothEarthHeights_diffractionModel(
            p,HTS_LIST[pathInd],HRS_LIST[pathInd]);
        
        const double LSPH = DiffractionLoss::se_diffLoss(p.back().d_km-p.front().d_km, 
            HTS_LIST[pathInd]+p.front().h_masl-height_eff_m.tx_val, //effective height relative to ground
            HRS_LIST[pathInd]+p.back().h_masl-height_eff_m.rx_val,
            ae, FREQ_MHZ_LIST[pathInd]/1000.0,
            0,Enumerations::PolarizationType::HorizontalPolarized); //Assume land, horizontal pol

        //suggested tolerance 0.1 dB from Notes page of spreadsheet. 
        //Validation data uses 3e8 for speed of light, model uses 2.998e8
        EXPECT_NEAR(EXPECTED_LSPH[pathInd],LSPH,0.1);
    }
}

TEST(DiffractionLossTests, DiffractionLoss_deltaBullingtonTest){
    // Arrange
	const std::vector<int> PATH_LIST = {
		1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4
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

    const std::vector<double> EXPECTED_LOSS = {
        34.44302652,42.32928281,36.34606901,27.69744396,42.99074274,
        69.83036947,90.68232743,66.46062522,46.38832061,51.38547874,
        17.62040631,25.1016907,20.25488733,10.1899085,20.73949248,
        7.870985135,0,15.30876397,6.526730443,0
    };

    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
        //input effective antanna heights relative to ground, frequency in GHz
        const PathProfile::Path p = PROFILE_LIST[PATH_LIST[pathInd]-1];
        const EffectiveEarth::HeightPair height_eff_m = EffectiveEarth::smoothEarthHeights_diffractionModel(
            p,HTS_LIST[pathInd],HRS_LIST[pathInd]);
        
        const double LOSS_VAL = DiffractionLoss::delta_bullington(
            p, 
            HTS_LIST[pathInd]+p.front().h_masl,
            HRS_LIST[pathInd]+p.back().h_masl,
            height_eff_m.tx_val,
            height_eff_m.rx_val,
            ae, FREQ_MHZ_LIST[pathInd]/1000.0,
            0,Enumerations::PolarizationType::HorizontalPolarized); //Assume land, horizontal pol

        //suggested tolerance 0.1 dB from Notes page of spreadsheet. 
        //Validation data uses 3e8 for speed of light, model uses 2.998e8
        EXPECT_NEAR(EXPECTED_LOSS[pathInd],LOSS_VAL,0.1);
    }
}