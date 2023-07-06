#include "DataGridHelpers.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

std::pair<MathHelpers::NeighborIntegerPair, MathHelpers::NeighborIntegerPair> DataGridHelpers::calculateBoundingBoxIntegerPairs(
			const GeodeticCoord& location, const double& gridResolution_deg, 
			const double& startLatitude_deg, const double& endLatitude_deg, 
			const double& startLongitude_deg, const double &endLongitude_deg) {

	double wrappedLongitude_deg = location.m_longitude_deg;
	// Edge case: When data grids are defined on a longitude axis of (0,360), instead of (-180,180)
	// Wrap GeodeticCoord object (which is defined from -180->180) so that it is defined from 0->360
	if (startLongitude_deg == 0.0 && endLongitude_deg == 360.0) {
		wrappedLongitude_deg = MathHelpers::unwrapValueAroundAxis(location.m_longitude_deg, 0.0, 360.0);
	}
	
	// NOTE: Start and end latitudes may be in descending order, if the DataGrid values are organized in descending order
	const double MIN_LATITUDE_DEG = std::min(startLatitude_deg, endLatitude_deg);
	const double MAX_LATITUDE_DEG = std::max(startLatitude_deg, endLatitude_deg);

	if (location.m_latitude_deg > MAX_LATITUDE_DEG || location.m_latitude_deg < MIN_LATITUDE_DEG) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: DataGridHelpers::calculateBoundingBoxIntegerPairs(): " 
					<< "Given coordinate's latitude falls outside of valid bounds [" 
					<< MIN_LATITUDE_DEG << " deg, " << MAX_LATITUDE_DEG << " deg]: " 
					<< location.m_latitude_deg << " deg! ";
		throw std::runtime_error(oStrStream.str());
	}
	if (wrappedLongitude_deg > endLongitude_deg || wrappedLongitude_deg < startLongitude_deg) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: DataGridHelpers::calculateBoundingBoxIntegerPairs(): " << 
					"Given coordinate's longitude falls outside of valid bounds [" 
					<< startLongitude_deg << " deg, " << endLongitude_deg << " deg]: " 
					<< wrappedLongitude_deg << " deg! ";
		throw std::runtime_error(oStrStream.str());
	}

	// Longitude is defined in _data so that
	// * startLongitude degrees is on rowInd = 0
	// * endLongitude degrees is on last rowInd
	const double lonColExact = abs(wrappedLongitude_deg - startLongitude_deg) / gridResolution_deg;
	MathHelpers::NeighborIntegerPair lonColNeighborPair(lonColExact);
	// Latitude is arranged from startLatitude to endLatitude
	// * startLatitude degrees is on rowInd = 0
	// * endLatitude is on last rowInd
	const double latRowExact = abs(location.m_latitude_deg - startLatitude_deg) / gridResolution_deg;
	MathHelpers::NeighborIntegerPair latRowNeighborPair(latRowExact);

	const double HIGH_POINT_LON_DEG = lonColNeighborPair.highPoint * gridResolution_deg + startLongitude_deg;

	// Don't allow the high longitude value to equal 360 degrees. Instead, wrap it around to the column for 0 degrees. 
	if (HIGH_POINT_LON_DEG == 360.0) {
		lonColNeighborPair.highPoint = 0.0;
	}

	return std::make_pair(lonColNeighborPair, latRowNeighborPair);
}

std::vector<DataGridHelpers::BoundingBoxGridPoint> DataGridHelpers::calculateBoundingBoxGridPointList(const GeodeticCoord& location, 
            const std::pair<MathHelpers::NeighborIntegerPair, MathHelpers::NeighborIntegerPair> &boundingBoxNeighborIntegerPairs,
	        const double& gridResolution_deg, const double &startLatitude_deg, 
            const double &startLongitude_deg, const bool &isLatitudeAscending) {
	const auto LON_COL_NEIGHBOR_PAIR = boundingBoxNeighborIntegerPairs.first;
	const auto LAT_ROW_NEIGHBOR_PAIR = boundingBoxNeighborIntegerPairs.second;

	// Convert back to geodetic
	double lon0_deg = startLongitude_deg + LON_COL_NEIGHBOR_PAIR.lowPoint * gridResolution_deg;
	double lon1_deg = startLongitude_deg + LON_COL_NEIGHBOR_PAIR.highPoint * gridResolution_deg;
	// Wrap these longitudes so that they are defined from -180->180, instead of 0->360 (otherwise GeodeticCoord will throw an error)
	if (lon0_deg > 180.0) {
		lon0_deg = MathHelpers::unwrapValueAroundAxis(lon0_deg, -180.0, 180.0);
	}
	if (lon1_deg > 180.0) {
		lon1_deg = MathHelpers::unwrapValueAroundAxis(lon1_deg, -180.0, 180.0);
	}

	double latitudeScalar = 1.0;
	if (!isLatitudeAscending) {
		latitudeScalar = -1.0;
	}

	const double lat0_deg = startLatitude_deg + LAT_ROW_NEIGHBOR_PAIR.lowPoint * latitudeScalar * gridResolution_deg;
	const double lat1_deg = startLatitude_deg + LAT_ROW_NEIGHBOR_PAIR.highPoint * latitudeScalar * gridResolution_deg;
	// Calculate weights for each point
	const double lon1Weight = LON_COL_NEIGHBOR_PAIR.weightFactor;
	const double lon0Weight = 1.0 - lon1Weight;
	const double lat1Weight = LAT_ROW_NEIGHBOR_PAIR.weightFactor;
	const double lat0Weight = 1.0 - lat1Weight;

	// Store data
	std::vector<BoundingBoxGridPoint> boundingBoxList;
	boundingBoxList.push_back(BoundingBoxGridPoint(GeodeticCoord(lon0_deg, lat0_deg, location.height_km), lon0Weight * lat0Weight));
	boundingBoxList.push_back(BoundingBoxGridPoint(GeodeticCoord(lon1_deg, lat0_deg, location.height_km), lon1Weight * lat0Weight));
	boundingBoxList.push_back(BoundingBoxGridPoint(GeodeticCoord(lon0_deg, lat1_deg, location.height_km), lon0Weight * lat1Weight));
	boundingBoxList.push_back(BoundingBoxGridPoint(GeodeticCoord(lon1_deg, lat1_deg, location.height_km), lon1Weight * lat1Weight));

	return boundingBoxList;
}

/// <summary>
/// Perform a 2D interpolation between 4 bounding box points, given the values associated with each points and weights associated with each axis
/// </summary>
/// <param name="valueList">Values associated with each point</param>
/// <param name="rowWeight">Weight associated with the upper value of the row axis (latitude)</param>
/// <param name="columnWeight">Weight associated with the upper value of the column axis (longitude)</param>
/// <returns>Interpolated value between the 4 given values</returns>
double DataGridHelpers::interpolate2D(const std::vector<double>& valueList, 
            const double& rowWeight, const double& columnWeight) {
	if (valueList.size() != 4) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: DataGridHelpers::interpolate2D(): " 
					<< "The provided list of values must have 4 elements: " << valueList.size() << "!";
		throw std::domain_error(oStrStream.str());
	}

	const double value00 = valueList[0];
	const double value01 = valueList[1];
	const double value10 = valueList[2];
	const double value11 = valueList[3];

	const double INTERP_ROW0 = MathHelpers::interpolate1D(value00, value01, columnWeight);
	const double INTERP_ROW1 = MathHelpers::interpolate1D(value10, value11, columnWeight);

	const double INTERP_COL = MathHelpers::interpolate1D(INTERP_ROW0, INTERP_ROW1, rowWeight);

	return INTERP_COL;
}

/// <summary>
/// Perform a 2D interpolation between 4 bounding box points, given the values associated with each points and some weights to scale those values
/// </summary>
/// <param name="valueList">Values associated with each point</param>
/// <param name="weightList">Weights to scale each value based on some custom logic</param>
/// <param name="rowWeight">Weight associated with the upper value of the row axis (latitude)</param>
/// <param name="columnWeight">Weight associated with the upper value of the column axis (longitude)</param>
/// <returns>Interpolated value between the 4 given values</returns>
double DataGridHelpers::interpolate2D(const std::vector<double>& valueList, 
            const std::vector<double>& weightList, const double& rowWeight, const double& columnWeight) {
	if (valueList.size() != 4 || weightList.size() != 4) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: DataGridHelpers::interpolate2D(): " 
					<< "The provided lists of values and weights must have 4 elements each: Num values = " 
					<< valueList.size() << ", Num weights = " << weightList.size() << "!";
		throw std::domain_error(oStrStream.str());
	}

	const double value00 = valueList[0] * weightList[0];
	const double value01 = valueList[1] * weightList[1];
	const double value10 = valueList[2] * weightList[2];
	const double value11 = valueList[3] * weightList[3];

	const double INTERP_ROW0 = MathHelpers::interpolate1D(value00, value01, columnWeight);
	const double INTERP_ROW1 = MathHelpers::interpolate1D(value10, value11, columnWeight);

	const double INTERP_COL = MathHelpers::interpolate1D(INTERP_ROW0, INTERP_ROW1, rowWeight);

	return INTERP_COL;
}