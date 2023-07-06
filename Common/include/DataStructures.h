#ifndef ITU_DATA_STRUCTURES_H
#define ITU_DATA_STRUCTURES_H

#include "MathHelpers.h"
#include "Enumerations.h"
#include "GeodeticCoord.h"
#include "ExceedanceDataGrid.h"

#include <vector>

namespace DataStructures {
	/* Represent each row of the oxygen coefficient table. */
	struct OxygenAttenData {
		OxygenAttenData(const double &freq_GHz, 
								const double &val1, const double& val2, const double& val3, 
								const double& val4, const double& val5, const double& val6)
			: _freq_GHz(freq_GHz), _a1(val1), _a2(val2), _a3(val3), _a4(val4), _a5(val5), _a6(val6) {}

		/// Frequency ID for the associated coefficients from the table
		double _freq_GHz;
		/// From Table 1 "Spectroscopic Data for Oxygen Attenuation", pg.8, ITU-R P.676-12
		double _a1, _a2, _a3, _a4, _a5, _a6;
	};

	static const uint32_t OXYGEN_COEFFS_TABLE_SIZE = 44;
	/// Table 1 from ITU-R P.676-12
	static const std::vector<OxygenAttenData> OXYGEN_COEFFS_TABLE {
		OxygenAttenData(50.474214,	0.975,	9.651,	6.69,	0.0,	2.566,	6.850),
		OxygenAttenData(50.987745,	2.529,	8.653,	7.17,	0.0,	2.246,	6.800),
		OxygenAttenData(51.50336,	6.193,	7.709,	7.64,	0.0,	1.947,	6.729),
		OxygenAttenData(52.021429,	14.32,	6.819,	8.11,	0.0,	1.667,	6.640),
		OxygenAttenData(52.542418,	31.24,	5.983,	8.58,	0.0,	1.388,	6.526),
		OxygenAttenData(53.066934,	64.29,	5.201,	9.06,	0.0,	1.349,	6.206),
		OxygenAttenData(53.595775,	124.6,	4.474,	9.55,	0.0,	2.227,	5.085),
		OxygenAttenData(54.130025,	227.3,	3.8,	9.96,	0.0,	3.170,	3.750),
		OxygenAttenData(54.67118,	389.7,	3.182,	10.37,	0.0,	3.558,	2.654),
		OxygenAttenData(55.221384,	627.1,	2.618,	10.89,	0.0,	2.560,	2.952),
		OxygenAttenData(55.783815,	945.3,	2.109,	11.34,	0.0,	-1.172,	6.135),
		OxygenAttenData(56.264774,	543.4,	0.014,	17.03,	0.0,	3.525,	-0.978),
		OxygenAttenData(56.363399,	1331.8,	1.654,	11.89,	0.0,	-2.378,	6.547),
		OxygenAttenData(56.968211,	1746.6,	1.255,	12.23,	0.0,	-3.545,	6.451),
		OxygenAttenData(57.612486,	2120.1,	0.91,	12.62,	0.0,	-5.416,	6.056),
		OxygenAttenData(58.323877,	2363.7,	0.621,	12.95,	0.0,	-1.932,	0.436),
		OxygenAttenData(58.446588,	1442.1,	0.083,	14.91,	0.0,	6.768,	-1.273),
		OxygenAttenData(59.164204,	2379.9,	0.387,	13.53,	0.0,	-6.561,	2.309),
		OxygenAttenData(59.590983,	2090.7,	0.207,	14.08,	0.0,	6.957,	-0.776),
		OxygenAttenData(60.306056,	2103.4,	0.207,	14.15,	0.0,	-6.395,	0.699),
		OxygenAttenData(60.434778,	2438.0,	0.386,	13.39,	0.0,	6.342,	-2.825),
		OxygenAttenData(61.150562,	2479.5,	0.621,	12.92,	0.0,	1.014,	-0.584),
		OxygenAttenData(61.800158,	2275.9,	0.91,	12.63,	0.0,	5.014,	-6.619),
		OxygenAttenData(62.41122,	1915.4,	1.255,	12.17,	0.0,	3.029,	-6.759),
		OxygenAttenData(62.486253,	1503.0,	0.083,	15.13,	0.0,	-4.499,	0.844),
		OxygenAttenData(62.997984,	1490.2,	1.654,	11.74,	0.0,	1.856,	-6.675),
		OxygenAttenData(63.568526,	1078.0,	2.108,	11.34,	0.0,	0.658,	-6.139),
		OxygenAttenData(64.127775,	728.7,	2.617,	10.88,	0.0,	-3.036,	-2.895),
		OxygenAttenData(64.67891,	461.3,	3.181,	10.38,	0.0,	-3.968,	-2.590),
		OxygenAttenData(65.224078,	274.0,	3.8,	9.96,	0.0,	-3.528,	-3.680),
		OxygenAttenData(65.764779,	153.0,	4.473,	9.55,	0.0,	-2.548,	-5.002),
		OxygenAttenData(66.302096,	80.4,	5.2,	9.06,	0.0,	-1.660,	-6.091),
		OxygenAttenData(66.836834,	39.8,	5.982,	8.58,	0.0,	-1.680,	-6.393),
		OxygenAttenData(67.369601,	18.56,	6.818,	8.11,	0.0,	-1.956,	-6.475),
		OxygenAttenData(67.900868,	8.172,	7.708,	7.64,	0.0,	-2.216,	-6.545),
		OxygenAttenData(68.431006,	3.397,	8.652,	7.17,	0.0,	-2.492,	-6.600),
		OxygenAttenData(68.960312,	1.334,	9.65,	6.69,	0.0,	-2.773,	-6.650),
		OxygenAttenData(118.750334,	940.3,	0.01,	16.64,	0.0,	-0.439,	0.079),
		OxygenAttenData(368.498246,	67.4,	0.048,	16.4,	0.0,	0.000,	0.000),
		OxygenAttenData(424.76302,	637.7,	0.044,	16.4,	0.0,	0.000,	0.000),
		OxygenAttenData(487.249273,	237.4,	0.049,	16.0,	0.0,	0.000,	0.000),
		OxygenAttenData(715.392902,	98.1,	0.145,	16.0,	0.0,	0.000,	0.000),
		OxygenAttenData(773.83949,	572.3,	0.141,	16.2,	0.0,	0.000,	0.000),
		OxygenAttenData(834.145546,	183.1,	0.145,	14.7,	0.0,	0.000,	0.000)
	};

	/* Represent each row of the water vapor coefficient table. */
	struct WaterVaporAttenData {
		WaterVaporAttenData(const double& freq_GHz,
			const double& val1, const double& val2, const double& val3,
			const double& val4, const double& val5, const double& val6)
			: _freq_GHz(freq_GHz), _b1(val1), _b2(val2), _b3(val3), _b4(val4), _b5(val5), _b6(val6) {}

		/// Frequency ID for the associated coefficients from the table
		double _freq_GHz;
		/// From Table 2 "Spectroscopic Data for Water Vapour Attenuation", pg.9, ITU-R P.676-12
		double _b1, _b2, _b3, _b4, _b5, _b6;
	};

	static const uint32_t WATER_COEFFS_TABLE_SIZE = 35;
	/// Table 2 from ITU-R P.676-12
	static const std::vector<WaterVaporAttenData> WATER_COEFFS_TABLE {
		WaterVaporAttenData(22.23508,	0.1079,	2.144,	26.38,	0.76,	5.087,	1.0),
		WaterVaporAttenData(67.80396,	0.0011,	8.732,	28.58,	0.69,	4.93,	0.82),
		WaterVaporAttenData(119.99594,	0.0007,	8.353,	29.48,	0.7,	4.78,	0.79),
		WaterVaporAttenData(183.310087,	2.273,	0.668,	29.06,	0.77,	5.022,	0.85),
		WaterVaporAttenData(321.22563,	0.047,	6.179,	24.04,	0.67,	4.398,	0.54),
		WaterVaporAttenData(325.152888,	1.514,	1.541,	28.23,	0.64,	4.893,	0.74),
		WaterVaporAttenData(336.227764,	0.001,	9.825,	26.93,	0.69,	4.74,	0.61),
		WaterVaporAttenData(380.197353,	11.67,	1.048,	28.11,	0.54,	5.063,	0.89),
		WaterVaporAttenData(390.134508,	0.0045,	7.347,	21.52,	0.63,	4.81,	0.55),
		WaterVaporAttenData(437.346667,	0.0632,	5.048,	18.45,	0.6,	4.23,	0.48),
		WaterVaporAttenData(439.150807,	0.9098,	3.595,	20.07,	0.63,	4.483,	0.52),
		WaterVaporAttenData(443.018343,	0.192,	5.048,	15.55,	0.6,	5.083,	0.5),
		WaterVaporAttenData(448.001085,	10.41,	1.405,	25.64,	0.66,	5.028,	0.67),
		WaterVaporAttenData(470.888999,	0.3254,	3.597,	21.34,	0.66,	4.506,	0.65),
		WaterVaporAttenData(474.689092,	1.26,	2.379,	23.2,	0.65,	4.804,	0.64),
		WaterVaporAttenData(488.490108,	0.2529,	2.852,	25.86,	0.69,	5.201,	0.72),
		WaterVaporAttenData(503.568532,	0.0372,	6.731,	16.12,	0.61,	3.98,	0.43),
		WaterVaporAttenData(504.482692,	0.0124,	6.731,	16.12,	0.61,	4.01,	0.45),
		WaterVaporAttenData(547.67644,	0.9785,	0.158,	26.0,	0.7,	4.5,	1.0),
		WaterVaporAttenData(552.02096,	0.184,	0.158,	26.0,	0.7,	4.5,	1.0),
		WaterVaporAttenData(556.935985,	497.0,	0.159,	30.86,	0.69,	4.552,	1.0),
		WaterVaporAttenData(620.700807,	5.015,	2.391,	24.38,	0.71,	4.856,	0.68),
		WaterVaporAttenData(645.766085,	0.0067,	8.633,	18.0,	0.6,	4.0,	0.5),
		WaterVaporAttenData(658.00528,	0.2732,	7.816,	32.1,	0.69,	4.14,	1.0),
		WaterVaporAttenData(752.033113,	243.4,	0.396,	30.86,	0.68,	4.352,	0.84),
		WaterVaporAttenData(841.051732,	0.0134,	8.177,	15.9,	0.33,	5.76,	0.45),
		WaterVaporAttenData(859.965698,	0.1325,	8.055,	30.6,	0.68,	4.09,	0.84),
		WaterVaporAttenData(899.303175,	0.0547,	7.914,	29.85,	0.68,	4.53,	0.9),
		WaterVaporAttenData(902.611085,	0.0386,	8.429,	28.65,	0.7,	5.1,	0.95),
		WaterVaporAttenData(906.205957,	0.1836,	5.11,	24.08,	0.7,	4.7,	0.53),
		WaterVaporAttenData(916.171582,	8.4,	1.441,	26.73,	0.7,	5.15,	0.78),
		WaterVaporAttenData(923.112692,	0.0079,	10.293,	29.0,	0.7,	5.0,	0.8),
		WaterVaporAttenData(970.315022,	9.009,	1.919,	25.5,	0.64,	4.94,	0.67),
		WaterVaporAttenData(987.926764,	134.6,	0.257,	29.85,	0.68,	4.55,	0.9),
		WaterVaporAttenData(1780.0,	17506.0,	0.952,	196.3,	2.0,	24.15,	5.0) 
	};

	// Table 3 on page #25 of ITU-R P.676-12
	static const std::vector<double> OxygenConstList = { 0.1597, 0.1066, 0.1325, 0.1242, 0.0938, 0.1448, 0.1374 };
	static const std::vector<double> OxygenFreqList = { 118.750334, 368.498246, 424.763020, 487.249273, 715.392902, 773.839490, 834.145546 };

	// Table 4 on page #26 of ITU-R P.676-12
	static const std::vector<double> WaterFreqList = { 22.23508, 183.310087, 325.152888, 380.197353, 439.150807, 448.001085, 474.689092,
	488.490108, 556.935985, 620.70087, 752.033113, 916.171582, 970.315022, 987.926764 };
	static const std::vector<double> WaterAList = { 1.52, 7.62, 1.56, 4.15, 0.20, 1.63, 0.76,
	0.26, 7.81, 1.25, 16.2, 1.47, 1.36, 1.60 };
	static const std::vector<double> WaterBList = { 2.56, 10.2, 2.70, 5.70, 0.91, 2.46, 2.22,
	2.49, 10.0, 2.35, 20.0, 2.58, 2.44, 1.86 };

	/// Days in each month for using monthly maps (ITU-R P.837-7)
	static const std::vector<double> DAYS_PER_MONTH = { 31,28.25,31,30,31,30,31,31,30,31,30,31 };

	/// Path loss scale factor in units of 1/(GHz * m) from 4 * PI / C
	static const double PATH_LOSS_SCALE_FACTOR = 41916.9;

	/// Margin for unit-fraction values
	static const double UNIT_FRACTION_MARGIN = 1.0e-5;

	/// Minimum allowed atmospheric (cloud and gas) exceedance is 1% (any exceedance lower than this value is accounted for by rain exceedance)
	static const double ATMOSPHERIC_EXCEEDANCE_MIN = 0.01;

	/// List of standard ITU probability data layers
	static const std::vector<ExceedanceDataGrid::ExceedanceLayer> ITU_LAYER_LIST {
		ExceedanceDataGrid::ExceedanceLayer("01", 0.1),
		ExceedanceDataGrid::ExceedanceLayer("02", 0.2),
		ExceedanceDataGrid::ExceedanceLayer("03", 0.3),
		ExceedanceDataGrid::ExceedanceLayer("05", 0.5),
		ExceedanceDataGrid::ExceedanceLayer("1", 1.),
		ExceedanceDataGrid::ExceedanceLayer("2", 2.),
		ExceedanceDataGrid::ExceedanceLayer("3", 3.),
		ExceedanceDataGrid::ExceedanceLayer("5", 5.),
		ExceedanceDataGrid::ExceedanceLayer("10", 10.),
		ExceedanceDataGrid::ExceedanceLayer("20", 20.),
		ExceedanceDataGrid::ExceedanceLayer("30", 30.),
		ExceedanceDataGrid::ExceedanceLayer("50", 50.),
		ExceedanceDataGrid::ExceedanceLayer("60", 60.),
		ExceedanceDataGrid::ExceedanceLayer("70", 70.),
		ExceedanceDataGrid::ExceedanceLayer("80", 80.),
		ExceedanceDataGrid::ExceedanceLayer("90", 90.),
		ExceedanceDataGrid::ExceedanceLayer("95", 95.),
		ExceedanceDataGrid::ExceedanceLayer("99", 99.)
	};
} // end namespace DataStructures
#endif /* ITU_DATA_STRUCTURES_H */
