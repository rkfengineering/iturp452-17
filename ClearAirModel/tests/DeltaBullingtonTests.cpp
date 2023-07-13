#include "gtest/gtest.h"
#include "DiffractionLoss.h"

#include <filesystem>

//Validation data from ITU validation spreadsheet titled "delB_valid_temp.xlsx", page "outputs"
//embedded in ITU validation document titled "Validation Examples for the delta Bullington diffraction prediction method"

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

    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
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

//Not sure what polarization type to use or what fraction over sea
//not specified by validation data
/*
TEST(DiffractionLossTests, DiffractionLoss_sphericalEarthFTTest){
    // Arrange
	const std::vector<int> PATH_LIST = {
		1,1,1,1,1,2,2,2,2,2,3,
	};
    const std::vector<double> HTS_LIST = {
        30,50,20,40,70,
        30,50,20,40,70,
        70
    };
    const std::vector<double> HRS_LIST = {
        30,10,20,50,5,
        30,10,20,50,5,
        5
    };
    const std::vector<double> FREQ_MHZ_LIST = {
        1000,2500,600,200,150,
        1000,2500,600,200,150,
        150
    };

    const std::vector<double> EXPECTED_FX = {
        -59.40894528,-85.5189285,-48.12283666,-34.24816912,-26.15770548,
        -102.9062648,-145.2758435,-84.49345195,-54.44670331,-48.32185499,
        -30.1070985
    };
    const std::vector<double> EXPECTED_GY1 = {
        39.52373705,62.46549974,30.34301718,18.03764354,16.59625588,
        7.40810458,23.83868156,0.127517272,-1.774893492,1.161982247,
        22.3511602
    };
    const std::vector<double> EXPECTED_GY2 = {
        4.536624536,-0.643878137,-3.134029165,-0.872488481,-23.57738733,
        26.04381871,35.89852443,18.04945817,10.466862,3.806254354,
        -15.23734738
    };
    
    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
        //hts,hrs needs to be in m asl, frequency in GHz
        const ProfilePath p = PROFILE_LIST[PATH_LIST[pathInd]-1];
        const double FT = DiffractionLoss::se_first_term(p.d.back()-p.d.front(), 
            HTS_LIST[pathInd]+p.h.front(),
            HRS_LIST[pathInd]+p.h.back(),
            ae, FREQ_MHZ_LIST[pathInd]/1000.0,
            0,Enumerations::PolarizationType::HorizontalPolarized); //Assume land, horizontal pol
        const double EXPECTED_FT = -EXPECTED_FX[pathInd] - EXPECTED_GY1[pathInd] - EXPECTED_GY2[pathInd];

        //suggested tolerance 0.1 dB from Notes page of spreadsheet. 
        //Validation data uses 3e8 for speed of light, model uses 2.998e8
        EXPECT_NEAR(EXPECTED_FT,FT,0.1);
    }
}
*/
//TEST(DiffractionLossTests, DiffractionLoss_sphericalEarthTest){
//TEST(DiffractionLossTests, DiffractionLoss_deltaBullingtonTest){