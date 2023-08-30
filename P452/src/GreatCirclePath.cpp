#include <cmath>
#include "P452/GreatCirclePath.h"
#include "Common/MathHelpers.h"

namespace{
    constexpr double kEarthRadius_m = 6371.0e3; // Volumetric mean radius of Earth
}

P452::GreatCirclePath::GreatCirclePath(const LatLonCoord &startLatLon, const LatLonCoord &endLatLon) {
	// Store LatLons
	startLatLon_ = startLatLon;
	endLatLon_ = endLatLon;

    // Convert LatLons into ECEFs
	startEcef_ = VectorHelpers::ConvertLatLonToECEF(startLatLon);
	endEcef_ = VectorHelpers::ConvertLatLonToECEF(endLatLon);

	// Calculate magnitude of start & end vectors for interpolation of path point heights later 
    startEcefMagnitude_m_ = VectorHelpers::CalculateMagnitudeMeters(startEcef_);
    endEcefMagnitude_m_ = VectorHelpers::CalculateMagnitudeMeters(endEcef_);

	totalAngle_rad_ = VectorHelpers::CalculateAngleBetweenRadians(startEcef_, endEcef_);
	VectorHelpers::ECEFCoordinate normalCrossVector = VectorHelpers::Normalize(
                VectorHelpers::CalculateCrossProduct(startEcef_,endEcef_)
    );

	//arbitrary perpendicular for special cases (edge case)
	if (totalAngle_rad_ == 0.0) {
		normalCrossVector = VectorHelpers::Normalize(VectorHelpers::ECEFCoordinate(-startEcef_.y_meters_, startEcef_.x_meters_, 0.0));
		//choose some path if points are on opposite sides of earth
		if(VectorHelpers::CalculateDotProduct(startEcef_, endEcef_) < 0.0) {
			totalAngle_rad_ = M_PI;
		}
	}
	normalStartVector_ = VectorHelpers::Normalize(startEcef_);
    normalCrossStartVector_ = VectorHelpers::Normalize(VectorHelpers::CalculateCrossProduct(normalCrossVector, normalStartVector_));
}

/// @brief Calculates a point at a given fraction along a known great circle path (using vector math)
/// @param fraction Fraction along the known great circle path (0-1)
/// @return Lat/lon point (accounts for height of start & end point) -- assumes spherical Earth
LatLonCoord P452::GreatCirclePath::calcPointAtFractionOfGreatCirclePath_vector(const double& fraction) const{

	double currentAngle_rad = totalAngle_rad_ * fraction;
	VectorHelpers::ECEFCoordinate resultVector = normalStartVector_ * cos(currentAngle_rad) + normalCrossStartVector_ * sin(currentAngle_rad);

	resultVector = resultVector * ItuModels::MathHelpers::interpolate1D(startEcefMagnitude_m_, endEcefMagnitude_m_, fraction);
	
	return VectorHelpers::ConvertECEFToLatLon(resultVector);
}

std::vector<LatLonCoord> P452::GreatCirclePath::calcPointsOnGreatCirclePath_vector(const uint32_t& numPoints) const{
	std::vector<LatLonCoord> pointsAlongPathList;

    const double fraction_denom = static_cast<double>(numPoints - 1);
    double fraction = 0.0;
	for (uint32_t pointInd = 0; pointInd < numPoints; pointInd++) {
        fraction = static_cast<double>(pointInd) / fraction_denom;
		pointsAlongPathList.push_back(calcPointAtFractionOfGreatCirclePath_vector(fraction));
	}
	return pointsAlongPathList;
}

std::vector<LatLonCoord> P452::GreatCirclePath::calcPointsOnGreatCirclePath_sphere(const uint32_t& numPoints) const{
    std::vector<LatLonCoord> pointsAlongPathList;
    
    const double startLat_rad = startLatLon_.lat_deg_ * 2.0 * M_PI / 180.0;
    const double startLon_rad = startLatLon_.lon_deg_ * 2.0 * M_PI / 180.0;
    const double endLat_rad = endLatLon_.lat_deg_ * 2.0 * M_PI / 180.0;
    const double endLon_rad = endLatLon_.lon_deg_ * 2.0 * M_PI / 180.0;

    const double deltaLat_rad = abs(startLat_rad - endLat_rad);
    const double deltaLon_rad = abs(startLon_rad - endLon_rad);
    const double distance_rad = sqrt(deltaLat_rad * deltaLat_rad + deltaLon_rad * deltaLon_rad);

    const double fraction_denom = static_cast<double>(numPoints - 1);
    double fraction = 0.0;
	for (uint32_t pointInd = 0; pointInd < numPoints; pointInd++) {
        fraction = static_cast<double>(pointInd) / fraction_denom;

        const double A = sin((1.0 - fraction) * distance_rad)/sin(distance_rad);
        const double B = sin(fraction * distance_rad) / sin(distance_rad);
        const double x = A*cos(startLat_rad)*cos(startLon_rad) + B*cos(endLat_rad)*cos(endLon_rad);
        const double y = A*cos(startLat_rad)*sin(startLon_rad) + B*cos(endLat_rad)*sin(endLon_rad);
        const double z = A*sin(startLat_rad) + B*sin(endLat_rad);
        
        const double interLat_rad = atan2(z,sqrt(x*x + y*y));
        const double interLon_rad = atan2(y,x);

        const double interLat_deg = interLat_rad * 180.0 / (2.0 * M_PI);
        const double interLon_deg = interLon_rad * 180.0 / (2.0 * M_PI);
        pointsAlongPathList.push_back(LatLonCoord(interLat_deg, interLon_deg));
    }

    return pointsAlongPathList;
}
