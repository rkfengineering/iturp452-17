#ifndef DATA_GRID_HELPERS_H
#define DATA_GRID_HELPERS_H

#include "GeodeticCoord.h"
#include "MathHelpers.h"

namespace DataGridHelpers {

struct BoundingBoxGridPoint {
	BoundingBoxGridPoint(const GeodeticCoord& location, const double& weight)
		: m_location(location), m_weight(weight) {}

	GeodeticCoord m_location;
	double m_weight;
};

std::pair<MathHelpers::NeighborIntegerPair, MathHelpers::NeighborIntegerPair> calculateBoundingBoxIntegerPairs(
			const GeodeticCoord& location, const double& gridResolution_deg, 
			const double& startLatitude_deg, const double& endLatitude_deg, 
			const double& startLongitude_deg, const double &endLongitude_deg);

std::vector<BoundingBoxGridPoint> calculateBoundingBoxGridPointList(const GeodeticCoord& location, 
			const std::pair<MathHelpers::NeighborIntegerPair, MathHelpers::NeighborIntegerPair> &boundingBoxNeighborIntegerPairs,
			const double& gridResolution_deg, const double &startLatitude_deg, 
			const double &startLongitude_deg, const bool &isLatitudeAscending = true);

/// <summary>
/// Perform a 2D interpolation between 4 bounding box points, given the values associated with each points and weights associated with each axis
/// </summary>
/// <param name="valueList">Values associated with each point</param>
/// <param name="rowWeight">Weight associated with the upper value of the row axis (latitude)</param>
/// <param name="columnWeight">Weight associated with the upper value of the column axis (longitude)</param>
/// <returns>Interpolated value between the 4 given values</returns>
double interpolate2D(const std::vector<double>& valueList, 
			const double& rowWeight, const double& columnWeight);

/// <summary>
/// Perform a 2D interpolation between 4 bounding box points, given the values associated with each points and some weights to scale those values
/// </summary>
/// <param name="valueList">Values associated with each point</param>
/// <param name="weightList">Weights to scale each value based on some custom logic</param>
/// <param name="rowWeight">Weight associated with the upper value of the row axis (latitude)</param>
/// <param name="columnWeight">Weight associated with the upper value of the column axis (longitude)</param>
/// <returns>Interpolated value between the 4 given values</returns>
double interpolate2D(const std::vector<double>& valueList, const std::vector<double>& weightList, 
			const double& rowWeight, const double& columnWeight);

} // end namespace

#endif /* DATA_GRID_HELPERS_H */