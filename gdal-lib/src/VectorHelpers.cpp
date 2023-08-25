#include "gdal-lib/VectorHelpers.h"

// Constants pulled from https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
namespace {
constexpr double kEquatorialRadius_meters = 6378137.0;
constexpr double kPolarRadius_meters = 6356752.3142;
} // Anonymous namespace


namespace VectorHelpers {

/*
* Converts a WGS-84 ECEF coordinate into a WGS-84 geodetic (lat/lon/height) coordinate using Ferrari's solution 
* (see: https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#The_application_of_Ferrari's_solution)
*/
LatLonCoord ConvertECEFToLatLon(const ECEFCoordinate& ecef) {
    const auto xMeters = ecef.x_meters_;
    const auto yMeters = ecef.y_meters_;
    const auto zMeters = ecef.z_meters_;

    // Initial setup
    const auto zSquared = zMeters * zMeters;
    const auto aSquared = kEquatorialRadius_meters * kEquatorialRadius_meters;
    const auto bSquared = kPolarRadius_meters * kPolarRadius_meters;
    const auto eSquared = (aSquared - bSquared) / aSquared;
    const auto ePrimeSquared = (aSquared - bSquared) / bSquared;

    // Calculate intermediate variables
    const auto p = sqrt(xMeters * xMeters + yMeters * yMeters);
    const auto F = 54.0 * bSquared * zSquared;
    const auto G = p * p + (1.0 - eSquared) * zSquared - eSquared * (aSquared - bSquared);
    const auto c = eSquared * eSquared * F * p * p / (G * G * G);
    const auto sTerm = 1.0 + c + sqrt(c * c + 2.0 * c);
    // Have to use C++ pow here because UNITS pow only supports int exponents
    const auto s = pow(sTerm, 1.0 / 3.0);
    const auto k = s + 1.0 + 1.0 / s;
    const auto P = F / (3.0 * k * k * G * G);
    const auto Q = sqrt(1.0 + 2.0 * eSquared * eSquared * P);
    const auto r0_term1 = -P * eSquared * p / (1.0 + Q);
    const auto r0_subTerm1 = 0.5 * aSquared * (1.0 + 1.0 / Q);
    const auto r0_subTerm2 = P * (1.0 - eSquared) * zSquared / (Q * (1.0 + Q));
    const auto r0_subTerm3 = 0.5 * P * p * p;
    const auto r0_term2 = sqrt(r0_subTerm1 - r0_subTerm2 - r0_subTerm3);
    const auto r0 = r0_term1 + r0_term2;
    const auto U_term1 = p - eSquared * r0;
    const auto U = sqrt(U_term1 * U_term1 + zSquared);
    const auto V = sqrt(U_term1 * U_term1 + (1.0 - eSquared) * zSquared);
    const auto z0 = bSquared * zMeters / (kEquatorialRadius_meters * V);

    // Calculate lat, lon, height
    const auto heightMeters = U * (1.0 - bSquared / (kEquatorialRadius_meters * V));
    const auto latRad = atan((zMeters + ePrimeSquared * z0) / p);
    const auto lonRad = atan2(yMeters, xMeters);

    const auto latDeg = latRad * 180.0 / M_PI;
    const auto lonDeg = lonRad * 180.0 / M_PI;

    return LatLonCoord(latDeg, lonDeg, heightMeters);
}

/*
* Converts a coordinate specified as WGS-84 geodetic (lat/lon/height) into a coordinate specified as WGS-84 Earth-Centered Earth-Fixed (ECEF)
*/
ECEFCoordinate ConvertLatLonToECEF(const LatLonCoord& latLon) {
    const auto latRad = latLon.lat_deg_ * M_PI / 180.0;
    const auto lonRad = latLon.lon_deg_ * M_PI / 180.0;

    const auto sinLatitude = sin(latRad);
    const auto cosLatitude = cos(latRad);
    const auto sinLongitude = sin(lonRad);
    const auto cosLongitude = cos(lonRad);

    const auto bSquared_over_aSquared =
        (kPolarRadius_meters * kPolarRadius_meters) / (kEquatorialRadius_meters * kEquatorialRadius_meters);
    const auto eSquared = 1.0 - bSquared_over_aSquared;
    // Calculate the prime vertical radius of Earth curvature
    const auto kN = kEquatorialRadius_meters / sqrt(1.0 - eSquared * sinLatitude * sinLatitude);

    // Convert lat/lon/height coordinate into X/Y/Z coordinate
    const double xMeters = (kN + latLon.height_meters_) * cosLatitude * cosLongitude;
    const double yMeters = (kN + latLon.height_meters_) * cosLatitude * sinLongitude;
    const double zMeters =
        (bSquared_over_aSquared * kN + latLon.height_meters_) * sinLatitude;

    // Represent final ECEF coordinate
    return ECEFCoordinate(xMeters, yMeters, zMeters);
}

double CalculateDistanceMeters(const ECEFCoordinate& startEcef, const ECEFCoordinate& endEcef) {
    const auto xDiff = endEcef.x_meters_ - startEcef.x_meters_;
    const auto xDiffSquared = xDiff * xDiff;
    const auto yDiff = endEcef.y_meters_ - startEcef.y_meters_;
    const auto yDiffSquared = yDiff * yDiff;
    const auto zDiff = endEcef.z_meters_ - startEcef.z_meters_;
    const auto zDiffSquared = zDiff * zDiff;

    return sqrt(xDiffSquared + yDiffSquared + zDiffSquared);
}

ECEFCoordinate Normalize(const ECEFCoordinate& ecefVector){
    const double vecMag = CalculateMagnitudeMeters(ecefVector);
    return ECEFCoordinate(
                ecefVector.x_meters_ / vecMag,
                ecefVector.y_meters_ / vecMag,
                ecefVector.z_meters_ / vecMag);
}

ECEFCoordinate CalculateCrossProduct(const ECEFCoordinate& firstVector, const ECEFCoordinate& secondVector){
    return ECEFCoordinate(firstVector.y_meters_ * secondVector.z_meters_ - firstVector.z_meters_ * secondVector.y_meters_, 
                            firstVector.z_meters_ * secondVector.x_meters_ - firstVector.x_meters_ * secondVector.z_meters_, 
                            firstVector.x_meters_ * secondVector.y_meters_ - firstVector.y_meters_ * secondVector.x_meters_);
}

double CalculateDotProduct(const ECEFCoordinate& firstVector, const ECEFCoordinate& secondVector) {
    return firstVector.x_meters_ * secondVector.x_meters_ + firstVector.y_meters_ * secondVector.y_meters_
           + firstVector.z_meters_ * secondVector.z_meters_;
}

double CalculateElevationDegrees(const ECEFCoordinate& earthEcef, const ECEFCoordinate& spaceEcef) {
    const double earthToSpaceVector_x = spaceEcef.x_meters_ - earthEcef.x_meters_;
    const double earthToSpaceVector_y = spaceEcef.y_meters_ - earthEcef.y_meters_;
    const double earthToSpaceVector_z = spaceEcef.z_meters_ - earthEcef.z_meters_;

    const ECEFCoordinate earthToSpaceVector(earthToSpaceVector_x, earthToSpaceVector_y, earthToSpaceVector_z);

    const auto earthToSpaceDotProduct = CalculateDotProduct(earthToSpaceVector, earthEcef);

    const ECEFCoordinate originEcef(0.0, 0.0, 0.0);
    const auto earthMagnitude = CalculateDistanceMeters(originEcef, earthEcef);
    const auto earthToSpaceMagnitude = CalculateDistanceMeters(earthEcef, spaceEcef);

    // Calculate elevation angle
    const double elevationAngle_rad =
        asin(earthToSpaceDotProduct / (earthMagnitude * earthToSpaceMagnitude));
    // Convert rad --> deg
    const double elevationAngle_deg = elevationAngle_rad * 180.0 / M_PI;
    return elevationAngle_deg;
}

double CalculateMagnitudeMeters(const ECEFCoordinate& ecefVector){
    return std::sqrt(ecefVector.x_meters_ * ecefVector.x_meters_ + ecefVector.y_meters_ * ecefVector.y_meters_
           + ecefVector.z_meters_ * ecefVector.z_meters_);
}

double CalculateAngleBetweenRadians(const ECEFCoordinate& firstVector, const ECEFCoordinate& secondVector){
    const double numer = CalculateDotProduct(firstVector,secondVector);
    const double denom = CalculateMagnitudeMeters(firstVector) * CalculateMagnitudeMeters(secondVector);
    return std::acos(numer / denom);
}

} // WGS84_VectorHelpers Namespace