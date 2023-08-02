#include "gtest/gtest.h"
#include "DeltaBullingtonTests.h"
#include "ClearAirModel/DiffractionLoss.h"
#include "ClearAirModel/ClearAirModelHelpers.h"

//Validation data from ITU validation spreadsheet titled "delB_valid_temp.xlsx", page "outputs"
//embedded in ITU validation document titled "Validation Examples for the delta Bullington diffraction prediction method"
//https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/Validation%20examples%20for%20the%20delta%20Bullington%20diffraction%20prediction%20method.docx

namespace {
	// Use when expected an exact match
	double constexpr TOLERANCE = 1.0e-6;
}
//allocate memory for the list of paths
std::vector<PathProfile::Path> DeltaBullingtonTests::m_profile_list={};

namespace ClearAirModel{

TEST_F(DeltaBullingtonTests, calculateDiffractionModelSmoothEarthHeightsTest){
    // Arrange

	const std::vector<int> PATH_LIST = {
		1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4
	};
    // Height of Tx Antenna above ground level (m)
    const std::vector<double> HTS_LIST = {
        30,50,20,40,70,
        30,50,20,40,70,
        30,50,20,40,70,
        30,50,20,40,70
    };
    // Height of Rx Antenna above ground level (m)
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
    
    // Height of Tx Antenna above ground level in smooth earth model (m)
    const std::vector<double> EXPECTED_HTSP = {
        211.3366061,233.4945304,204.7165548,215.727222,252.3138569,37.97366242,
        57.97664284,27.97366242,47.97366242,78.00115495,110.1894677,132.6734925,
        102.8893122,118.7202148,152.2011247,139.8495825,159.8495825,129.8495825,149.8495825,179.8495825
    };
    // Height of Rx Antenna above ground level in smooth earth model (m)
    const std::vector<double> EXPECTED_HRSP = {
        30,10,20,50,5,119.2363492,99.82603022,109.2363492,139.2363492,99.67579634,30,
        11.9416848,22.52524312,50,5.664441408,543.1729418,523.1729418,533.1729418,
        563.1729418,518.1729418,
    };

    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
        const PathProfile::Path path = m_profile_list[PATH_LIST[pathInd]-1];
        //Calculate actual antenna heights relative to sea level
        const double HTS_AMSL = HTS_LIST[pathInd]+path.front().h_asl_m;
        const double HRS_AMSL = HRS_LIST[pathInd]+path.back().h_asl_m;

        const auto DiffractionModel = DiffractionLoss(
            path, 
            HTS_AMSL,
            HRS_AMSL,
            FREQ_MHZ_LIST[pathInd]/1000.0,
            53,//delta N isn't tested here
            Enumerations::PolarizationType::HorizontalPolarized,
            0.1, //p_percent isn't tested here
            3.0, //b0_percent isn't tested here
            path.calcFracOverSea()
        );

        //Calculate ground level height relative to sea level at tx,rx in smooth earth model 
        const auto [EFF_HEIGHT_TX, EFF_HEIGHT_RX] = DiffractionModel.calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m();

        //Equation 38 to calculate modified antenna heights
        EXPECT_NEAR(EXPECTED_HTSP[pathInd],HTS_AMSL-EFF_HEIGHT_TX,TOLERANCE);
        EXPECT_NEAR(EXPECTED_HRSP[pathInd],HRS_AMSL-EFF_HEIGHT_RX,TOLERANCE);
    }
}

TEST_F(DeltaBullingtonTests, DiffractionLoss_calcBullingtonLossTest){
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
    
    //Bullington loss of actual profile path
    const std::vector<double> EXPECTED_LBA = {
        32.84284766, 37.229036, 31.22595813, 24.87792613, 24.78223928,
        36.20419305,44.81240928,34.43917939,28.37368149,34.61746032,
        17.62040631,25.1016907,20.25488733,10.1899085,16.27593106,
        7.8709851352,0,15.30876397,6.526730443,0
    };

    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); ++pathInd) {
        //input heights need to be in m asl, frequency in GHz
        const PathProfile::Path path = m_profile_list[PATH_LIST[pathInd]-1];

        const auto DiffractionModel = DiffractionLoss(
            path, 
            HTS_LIST[pathInd] + path.front().h_asl_m,
            HRS_LIST[pathInd] + path.back().h_asl_m,
            FREQ_MHZ_LIST[pathInd]/1000.0,
            53,//delta N isn't tested here
            Enumerations::PolarizationType::HorizontalPolarized,
            0.1, //p_percent isn't tested here
            3.0, //b0_percent isn't tested here
            path.calcFracOverSea()
        );

        const double LBA = DiffractionModel.calcBullingtonLoss_dB(
            path,
            HTS_LIST[pathInd] + path.front().h_asl_m,
            HRS_LIST[pathInd] + path.back().h_asl_m,
            m_effEarthRadius_km
        );
        
        //suggested tolerance 0.1 dB from Notes page of spreadsheet. 
        //Validation data uses 3e8 for speed of light, model uses 2.998e8
        EXPECT_NEAR(EXPECTED_LBA[pathInd],LBA,m_deltaBullingtonTolerance_dB);
    }
}

//Two tests are excluded (row 6 and row 17) since they take a different branch in the code
//for the spherical earth calculations, resulting in a different modified effective earth radius being used
//for the intermediate calculations that produce the values used to validate this first term function
TEST_F(DeltaBullingtonTests, DiffractionLoss_calcSphericalEarthLossFirstTermTest){
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

    //Intermediate Values

    //Path distance effect on spherical diffraction loss
    const std::vector<double> EXPECTED_FX = {
        -59.40894528,-85.5189285,-48.12283666,-26.15770548,
        -102.9062648,-145.2758435,-84.49345195,-54.44670331,-48.32185499
    };
    //Tx height effect on spherical diffraction loss
    const std::vector<double> EXPECTED_GY1 = {
        39.52373705,62.46549974,30.34301718,16.59625588,
        7.40810458,23.83868156,0.127517272,-1.774893492,1.161982247
    };
    //Rx height effect on spherical diffraction loss
    const std::vector<double> EXPECTED_GY2 = {
        4.536624536,-0.643878137,-3.134029165,-23.57738733,
        26.04381871,35.89852443,18.04945817,10.466862,3.806254354
    };
    
    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
        //input effective antanna heights relative to ground, frequency in GHz
        const PathProfile::Path path = m_profile_list[PATH_LIST[pathInd]-1];

        const auto DiffractionModel = DiffractionLoss(
            path, 
            HTS_LIST[pathInd] + path.front().h_asl_m,
            HRS_LIST[pathInd] + path.back().h_asl_m,
            FREQ_MHZ_LIST[pathInd]/1000.0,
            53,//delta N isn't tested here
            Enumerations::PolarizationType::HorizontalPolarized,
            0.1, //p_percent isn't tested here
            3.0, //b0_percent isn't tested here
            path.calcFracOverSea()
        );
        
        //First term of spherical diffraction loss
        const double FT = DiffractionModel.calcSphericalEarthDiffraction_firstTerm_dB(m_effEarthRadius_km);

        //Expected First term is a sum of intermediate values
        const double EXPECTED_FT = -EXPECTED_FX[pathInd] - EXPECTED_GY1[pathInd] - EXPECTED_GY2[pathInd];

        EXPECT_NEAR(EXPECTED_FT,FT,m_deltaBullingtonTolerance_dB);
    }
}

TEST_F(DeltaBullingtonTests, DiffractionLoss_calcSphericalEarthLossTest){
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

    //Spherical Diffraction Loss (dB)
    const std::vector<double> EXPECTED_LSPH = {
        15.3485837,23.6973069,20.91384864,14.32933641,33.13883693,
        69.45434151,85.53863749,66.31647651,45.7547348,43.35361839,
        0,0,0,0,10.63765222,
        0,0,0,0,0
    };

    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
        //input effective antanna heights relative to ground, frequency in GHz
        const PathProfile::Path path = m_profile_list[PATH_LIST[pathInd]-1];

        const auto DiffractionModel = DiffractionLoss(
            path, 
            HTS_LIST[pathInd] + path.front().h_asl_m,
            HRS_LIST[pathInd] + path.back().h_asl_m,
            FREQ_MHZ_LIST[pathInd]/1000.0,
            53,//delta N isn't tested here
            Enumerations::PolarizationType::HorizontalPolarized,
            0.1, //p_percent isn't tested here
            3.0, //b0_percent isn't tested here
            path.calcFracOverSea()
        );

        const double LSPH = DiffractionModel.calcSphericalEarthDiffractionLoss_dB(m_effEarthRadius_km);
        
        EXPECT_NEAR(EXPECTED_LSPH[pathInd],LSPH,m_deltaBullingtonTolerance_dB);
    }
}

TEST_F(DeltaBullingtonTests, DiffractionLoss_calcDeltaBullingtonLossTest){
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

    //Loss predicted from the Delta-Bullington Model
    const std::vector<double> EXPECTED_LOSS = {
        34.44302652,42.32928281,36.34606901,27.69744396,42.99074274,
        69.83036947,90.68232743,66.46062522,46.38832061,51.38547874,
        17.62040631,25.1016907,20.25488733,10.1899085,20.73949248,
        7.870985135,0,15.30876397,6.526730443,0
    };

    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
        //input effective antanna heights relative to ground, frequency in GHz
        const PathProfile::Path path = m_profile_list[PATH_LIST[pathInd]-1];
        
        const auto DiffractionModel = DiffractionLoss(
            path, 
            HTS_LIST[pathInd] + path.front().h_asl_m,
            HRS_LIST[pathInd] + path.back().h_asl_m,
            FREQ_MHZ_LIST[pathInd]/1000.0,
            53,//delta N isn't tested here
            Enumerations::PolarizationType::HorizontalPolarized,
            0.1, //p_percent isn't tested here
            3.0, //b0_percent isn't tested here
            path.calcFracOverSea()
        );

        const double LOSS_VAL = DiffractionModel.calcDeltaBullingtonLoss_dB(m_effEarthRadius_km);
 
        //suggested tolerance 0.1 dB from Notes page of spreadsheet. 
        //Validation data uses 3e8 for speed of light, model uses 2.998e8
        EXPECT_NEAR(EXPECTED_LOSS[pathInd],LOSS_VAL,m_deltaBullingtonTolerance_dB);
    }
}

TEST_F(DeltaBullingtonTests, DiffractionLoss_calcSmoothEarthBullingtonLossTest){
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

    //Loss predicted from the Delta-Bullington Model
    const std::vector<double> EXPECTED_LOSS = {
        13.74840484,18.5970601,15.79373777,11.50981857,14.93033347,
        35.82816509,39.66871934,34.29503068,27.74009568,26.58559997,
        0,0,0,0,6.17409081,
        0,0,0,0,0
    };

    for (uint32_t pathInd = 0; pathInd < PATH_LIST.size(); pathInd++) {
        //input effective antanna heights relative to ground, frequency in GHz
        const PathProfile::Path path = m_profile_list[PATH_LIST[pathInd]-1];
        const double FREQ_GHZ = FREQ_MHZ_LIST[pathInd]/1000.0; //ERROR things break if we use this expression directly. probably cuz everything is done by reference
        const auto DiffractionModel = DiffractionLoss(
            path, 
            HTS_LIST[pathInd] + path.front().h_asl_m,
            HRS_LIST[pathInd] + path.back().h_asl_m,
            FREQ_GHZ,
            53,//delta N isn't tested here
            Enumerations::PolarizationType::HorizontalPolarized,
            0.1, //p_percent isn't tested here
            3.0, //b0_percent isn't tested here
            path.calcFracOverSea()
        );

        //modified heights and zero profile
        PathProfile::Path zeroHeightPath;
        zeroHeightPath.clear();
        for (auto point : path){
            zeroHeightPath.push_back(PathProfile::ProfilePoint(point.d_km, 0.0));
        }
        //effective heights for smooth path
        const auto [eff_height_itx_asl_m,eff_height_irx_asl_m] = 
                DiffractionModel.calcSmoothEarthTxRxHeights_DiffractionModel_amsl_m();
        const double mod_height_tx_asl_m= HTS_LIST[pathInd] + path.front().h_asl_m- eff_height_itx_asl_m;
        const double mod_height_rx_asl_m= HRS_LIST[pathInd] + path.back().h_asl_m- eff_height_irx_asl_m;

        const double LBS = DiffractionModel.calcBullingtonLoss_dB(
            zeroHeightPath,
            mod_height_tx_asl_m,
            mod_height_rx_asl_m,
            m_effEarthRadius_km
        );
 
        //suggested tolerance 0.1 dB from Notes page of spreadsheet. 
        //Validation data uses 3e8 for speed of light, model uses 2.998e8
        EXPECT_NEAR(EXPECTED_LOSS[pathInd],LBS,m_deltaBullingtonTolerance_dB);
    }
}

}//end namespace ClearAirModel
