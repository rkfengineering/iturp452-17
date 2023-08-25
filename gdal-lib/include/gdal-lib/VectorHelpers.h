#ifndef VECTOR_HELPERS_H
#define VECTOR_HELPERS_H

#include <cmath>
#include <iostream>

#include "LatLonCoord.h"

namespace VectorHelpers {

struct ECEFCoordinate {
        explicit ECEFCoordinate(const double& x_meters, const double& y_meters, const double& z_meters)
            : x_meters_(x_meters)
            , y_meters_(y_meters)
            , z_meters_(z_meters) {
        }

        ECEFCoordinate()
            : x_meters_{0.0}
            , y_meters_{0.0}
            , z_meters_{0.0} {
        }

        double x_meters_;
        double y_meters_;
        double z_meters_;
    
        //WARNING: NEW
        inline ECEFCoordinate operator * (const double scalar) const {
            return ECEFCoordinate(x_meters_ * scalar, y_meters_ * scalar, z_meters_ * scalar);		
        }
        inline ECEFCoordinate operator + (const ECEFCoordinate &other) const {
            return ECEFCoordinate(x_meters_ + other.x_meters_, y_meters_ + other.y_meters_,z_meters_ + other.z_meters_);
        }
        friend std::ostream& operator<< (std::ostream &out, ECEFCoordinate const& data) {
            out << "(x = " << data.x_meters_ << " m, y = " << data.y_meters_ << " m, z = " <<data.z_meters_ << " m)";
            return out;
        }
};

/*!
* 
* @param latLon A latitude, longitude, and height value struct to describe a WGS-84 coordinate
* @return Struct containing the WGS-84 coordinate in the XYZ ECEF coordinate frame
*/
ECEFCoordinate ConvertLatLonToECEF(const LatLonCoord& latLon);
/*!
* 
* @param ecefCoordinate An X, Y, and Z value struct to describe a WGS-84 coordinate
* @return Struct containing the WGS-84 coordinate in the LLA geodetic coordinate frame
*/
LatLonCoord ConvertECEFToLatLon(const ECEFCoordinate& ecef);

double CalculateDistanceMeters(const ECEFCoordinate& startEcef, const ECEFCoordinate& endEcef);

ECEFCoordinate Normalize(const ECEFCoordinate& ecefVector);
ECEFCoordinate CalculateCrossProduct(const ECEFCoordinate& firstVector, const ECEFCoordinate& secondVector);
double CalculateDotProduct(const ECEFCoordinate& firstVector, const ECEFCoordinate& secondVector);
double CalculateElevationDegrees(const ECEFCoordinate& earthEcef, const ECEFCoordinate& spaceEcef);

//WARNING: NEW
double CalculateMagnitudeMeters(const ECEFCoordinate& ecefVector);
double CalculateAngleBetweenRadians(const ECEFCoordinate& firstVector, const ECEFCoordinate& secondVector);

}// namespace WGS84_VectorHelpers

#endif// WGS84_VECTOR_HELPERS_H