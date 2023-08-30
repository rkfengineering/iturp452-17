#include "gtest/gtest.h"
#include "MainModel/PathProfile.h"
#include "MainModel/CalculationHelpers.h"
#include "MainModel/ClearAirModelHelpers.h"
#include "MainModel/DataGridTxt.h"
#include "MainModel/BasicProp.h"
#include "MainModel/TropoScatter.h"
#include "MainModel/DataLoader.h"

#include "Common/PowerUnitConversionHelpers.h"
#include "Common/DataStructures.h"
#include <filesystem>

//Example Profile Path from ITU validation spreadsheet titled "delB_valid_temp.xlsx", pages "Path 1" to "Path 4"
//embedded in ITU validation document titled "Validation Examples for the delta Bullington diffraction prediction method"

//Unit Tests for smaller functions (not directly using ITU Validation data)

namespace {
	// Use when expected an exact match
	double constexpr TOLERANCE = 1.0e-6;
	const std::filesystem::path clearAirPathsFullPath = CMAKE_CLEARAIR_SRC_DIR / std::filesystem::path("tests/test_paths");
	const std::filesystem::path clearAirDataFullPath = CMAKE_CLEARAIR_SRC_DIR / std::filesystem::path("data");
}

namespace ITUR_P452{

TEST(ProfilePathTests, loadProfileTest){
	const std::vector<std::string> PROFILE_LIST = {
		"dbull_path1.csv", "test_profile_mixed_109km.csv"
	};
	const std::vector<double> EXPECTED_FIRST_D_LIST = {
		0,0
	};
	const std::vector<double> EXPECTED_FIRST_H_LIST = {
		800,40
	};
	const std::vector<int> EXPECTED_FIRST_ZONE_LIST = {
		0,1
	};
	const std::vector<double> EXPECTED_LAST_D_LIST = {
		83.07993853,109
	};
	const std::vector<double> EXPECTED_LAST_H_LIST = {
		302,183
	};
	const std::vector<int> EXPECTED_LAST_ZONE_LIST = {
		0,2
	};
	const std::vector<double> EXPECTED_LENGTH_LIST = {
		1663,110
	};
	for (uint16_t profileInd = 0; profileInd < PROFILE_LIST.size(); profileInd++) {
		const PathProfile::Path PROFILE((clearAirPathsFullPath/std::filesystem::path(PROFILE_LIST[profileInd])).string());
		EXPECT_NEAR(EXPECTED_FIRST_D_LIST[profileInd], PROFILE.front().d_km, TOLERANCE);
		EXPECT_NEAR(EXPECTED_FIRST_H_LIST[profileInd], PROFILE.front().h_asl_m, TOLERANCE);
		EXPECT_EQ(EXPECTED_FIRST_ZONE_LIST[profileInd], static_cast<int>(PROFILE.front().zone));
		EXPECT_NEAR(EXPECTED_LAST_D_LIST[profileInd], PROFILE.back().d_km, TOLERANCE);
		EXPECT_NEAR(EXPECTED_LAST_H_LIST[profileInd], PROFILE.back().h_asl_m, TOLERANCE);
		EXPECT_EQ(EXPECTED_LAST_ZONE_LIST[profileInd], static_cast<int>(PROFILE.back().zone));
		EXPECT_NEAR(EXPECTED_LENGTH_LIST[profileInd], PROFILE.size(), TOLERANCE);
	}	
}

//Compared against NORM.INV Built-in Macro from Microsoft Excel
TEST(HelpersTests, InvCumNormTest){
	const std::vector<double> INPUTS_LIST = {
		1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5,
		0.01345, 0.42039, 0.0042598, 0.000050938
	};
	const std::vector<double> EXPECTED_RESULTS_LIST = {
		-4.753424309, -4.264890794, -3.719016485,
		-3.090232306, -2.326347874, -1.281551566, 0.0, -2.212965852,
		-0.200895866, -2.630752722, -3.886079871,
	};

	for (uint16_t ind = 0; ind < INPUTS_LIST.size(); ind++) {
		const double OUTPUT = CalculationHelpers::inv_cum_norm(INPUTS_LIST[ind]);
		//Maximum error of 0.00054 from Attachment 3 to Annex 1 ITU-R P.452-17
		EXPECT_NEAR(EXPECTED_RESULTS_LIST[ind], OUTPUT, 0.00054);
	}
	
	//use 1e-6 for all probability below 1e-6
	EXPECT_NEAR(CalculationHelpers::inv_cum_norm(1e-7), CalculationHelpers::inv_cum_norm(1e-6), TOLERANCE);
}

//Check frequency to wavelength conversion using values from 
//https://www.translatorscafe.com/unit-converter/en-US/frequency-wavelength/5-27/gigahertz-wavelength%20in%20metres/
TEST(HelpersTests, convertFreqGHzToWavelengthMTest){
	const std::vector<double> INPUT_FREQ_GHZ_LIST = {
		0.5,1.0,2.0,10.0
	};
	const std::vector<double> EXPECTED_WAVELENGTH_M_LIST = {
		0.599584916,0.299792458,0.149896229,0.0299792458
	};

	for (uint16_t ind = 0; ind < INPUT_FREQ_GHZ_LIST.size(); ind++) {
		const double OUTPUT = CalculationHelpers::convert_freqGHz_to_wavelength_m(INPUT_FREQ_GHZ_LIST[ind]);
		//Maximum error of 0.00054 from Attachment 3 to Annex 1 ITU-R P.452-17
		EXPECT_NEAR(EXPECTED_WAVELENGTH_M_LIST[ind], OUTPUT, TOLERANCE);
	}
}

//This helper function returns the endpoints of a least squares smooth earth approximation
//This is used in calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m to calculate the effective heights 
//for the smooth earth Bullington Loss after compensating for terrain obstructions and 
//setting a lower yheight limit at the actual terrain height 

//Endpoints estimated from Excel Linear Trendline
TEST(HelpersTests, calcLeastSquaresSmoothEarthHeightsHelper){
	const double EXPECTED_START = 635.2;
	const double EXPECTED_END = 357.3;

	const PathProfile::Path p((clearAirPathsFullPath/std::filesystem::path("dbull_path1.csv")).string());
	const auto [eff_height_tx, eff_height_rx] = Helpers::calcLeastSquaresSmoothEarthTxRxHeights_helper_amsl_m(p);

	//high tolerance since the values were read off a graph
	EXPECT_NEAR(EXPECTED_START,eff_height_tx,0.5);
	EXPECT_NEAR(EXPECTED_END,eff_height_rx,0.5);
}

//Compare free space path loss value against other existing implementation
//The constant used in Eq 8 uses less sig figs
TEST(BasicPropTests, calcFreeSpacePathLoss){
	const double INPUT_FREQ_GHz = 2;
	const double INPUT_DIST_KM = 500;

	const double EXPECTED_LOSS = ItuModels::PowerUnitConversionHelpers::convertWattsToDb(
		ItuModels::DataStructures::PATH_LOSS_SCALE_FACTOR * INPUT_DIST_KM * INPUT_FREQ_GHz) * 2.0;

	//arbitrary inputs to create the object
	const auto BasicPropModel = BasicProp(INPUT_DIST_KM, 0,0, INPUT_FREQ_GHz, 300, 1000, 
            0, 0.1, 3, ITUR_P452::TxRxPair{20,20});	
	
	const double VAL_LOSS =  BasicPropModel.calcFreeSpacePathLoss_dB(INPUT_DIST_KM,INPUT_FREQ_GHz);

	//Eq 8 Constant used has resolution up to 0.1 dB
	EXPECT_NEAR(EXPECTED_LOSS,VAL_LOSS,0.05);
}

//check to see if they load without error
TEST(DataGridTests, loadDataGridTest){
	const std::vector<std::string> DATA_LIST = {
		"DN50.TXT", "N050.TXT"
	};
	const double RESOLUTION = 1.5;
	const std::vector<double> EXPECTED_FIRST = {
		40.726,317.248
	};
	for (uint16_t dataInd = 0; dataInd < DATA_LIST.size(); dataInd++) {
		const DataGridTxt data((clearAirDataFullPath/std::filesystem::path(DATA_LIST[dataInd])).string(),RESOLUTION);

		EXPECT_NEAR(EXPECTED_FIRST[dataInd], data.interpolate2D(ItuModels::GeodeticCoord(0.0,89.999999)), TOLERANCE);
	}	
}

//Check values at lat/lon coordinates
//expected values checked by hand
TEST(DataGridTests, getDataGridValuesDN50Test){
	const double RESOLUTION = 1.5;
	const DataGridTxt data((clearAirDataFullPath/std::filesystem::path("DN50.TXT")).string(),RESOLUTION);

	const std::vector<double> INPUT_LON = {
		36,-4.5,179.99999999
	};
	const std::vector<double> INPUT_LAT = {
		61.5,24.0,-58.5
	};
	const std::vector<double> EXPECTED_DN50 = {
		36.530,40.412,38.788
	};

	for (uint16_t coordInd = 0; coordInd < INPUT_LON.size(); coordInd++) {
		EXPECT_NEAR(EXPECTED_DN50[coordInd], data.interpolate2D(ItuModels::GeodeticCoord(INPUT_LON[coordInd],INPUT_LAT[coordInd])), TOLERANCE);
		EXPECT_NEAR(EXPECTED_DN50[coordInd], DataLoader::fetchRadioRefractivityIndexLapseRate(ItuModels::GeodeticCoord(INPUT_LON[coordInd],INPUT_LAT[coordInd])), TOLERANCE);
	}	
}

TEST(DataGridTests, fetchDataGridValuesN050Test){
	const double RESOLUTION = 1.5;
	const DataGridTxt data((clearAirDataFullPath/std::filesystem::path("N050.TXT")).string(),RESOLUTION);

	const std::vector<double> INPUT_LON = {
		36,-4.5,179.99999999
	};
	const std::vector<double> INPUT_LAT = {
		61.5,24.0,-58.5
	};
	const std::vector<double> EXPECTED_N050 = {
		316.375,314.726,317.199
	};

	for (uint16_t coordInd = 0; coordInd < INPUT_LON.size(); coordInd++) {
		EXPECT_NEAR(EXPECTED_N050[coordInd], data.interpolate2D(ItuModels::GeodeticCoord(INPUT_LON[coordInd],INPUT_LAT[coordInd])), TOLERANCE);
		EXPECT_NEAR(EXPECTED_N050[coordInd], DataLoader::fetchSeaLevelSurfaceRefractivity(ItuModels::GeodeticCoord(INPUT_LON[coordInd],INPUT_LAT[coordInd])), TOLERANCE);
	}	
}

}//end namespace ITUR_P452