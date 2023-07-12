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

static const std::vector<ProfilePath> PROFILE_LIST = {
    ProfilePath(testPathFileDir/std::filesystem::path("path1.csv")),
    ProfilePath(testPathFileDir/std::filesystem::path("path2.csv")),
    ProfilePath(testPathFileDir/std::filesystem::path("path3.csv")),
    ProfilePath(testPathFileDir/std::filesystem::path("path4.csv")),
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
        const ProfilePath p = PROFILE_LIST[PATH_LIST[pathInd]-1];
        const double LBA = DiffractionLoss::bullLoss(p.d, p.h,
            HTS_LIST[pathInd]+p.h.front(),
            HRS_LIST[pathInd]+p.h.back(),
            ae, FREQ_MHZ_LIST[pathInd]/1000.0);

        //suggested tolerance 0.1 dB from Notes page of spreadsheet. 
        //Validation data uses 3e8 for speed of light, model uses 2.998e8
        EXPECT_NEAR(EXPECTED_LBA[pathInd],LBA,0.1);
    }
}

//Not sure what polarization type to use or what fraction over sea
//not specified by validation data
//TEST(DiffractionLossTests, DiffractionLoss_sphericalEarthFTTest){
//TEST(DiffractionLossTests, DiffractionLoss_sphericalEarthTest){
//TEST(DiffractionLossTests, DiffractionLoss_deltaBullingtonTest){