#ifndef GEODETIC_COORD_H
#define GEODETIC_COORD_H

#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>

/// <summary>
/// Contains a 3D Earth-fixed geodetic coordinate.
/// This is in the WGS-84 ellipsoid, so any conversion functions must follow the WGS-84 conventions.
/// NOTE: In some cases, the height is unused and therefore it has been made an optional parameter.
/// </summary>
class GeodeticCoord {
public:
    GeodeticCoord()
                : m_longitude_deg(std::numeric_limits<double>::quiet_NaN()), 
                m_latitude_deg(std::numeric_limits<double>::quiet_NaN()), 
                height_km(0.0) {}

    GeodeticCoord(const double& lon_deg, const double& lat_deg, double height_km = 0.0)
                : height_km(height_km) {
        // Enforce -90 --> 90 deg bounds for latitude
        if (lat_deg >= 90.0 || lat_deg <= -90.0) {
            std::ostringstream oStrStream;
            oStrStream << std::fixed;
            oStrStream << std::setprecision(4);
            oStrStream << "ERROR: GeodeticCoord::GeodeticCoord(): " 
                        << "Given coordinate's latitude falls outside of valid bounds [-90 deg, +90 deg]: " << lat_deg << " deg! ";
            throw std::runtime_error(oStrStream.str());
        }
        // Enforce -180 --> 180 deg bounds for longitude
        if (lon_deg >= 180.0 || lon_deg <= -180.0) {
            std::ostringstream oStrStream;
            oStrStream << std::fixed;
            oStrStream << std::setprecision(4);
            oStrStream << "ERROR: GeodeticCoord::GeodeticCoord(): " 
                        << "Given coordinate's longitude falls outside of valid bounds [-180 deg, +180 deg]: " << lon_deg << " deg! ";
            throw std::runtime_error(oStrStream.str());
        }

        m_latitude_deg = lat_deg;
        m_longitude_deg = lon_deg;
    }

    /// <summary>
    /// Check if the contents of this objects are NaN
    /// </summary>
    /// <returns>Boolean indicating whether any contents of this object are NaN</returns>
    bool isNull() const;

    /// Longitude referenced to WGS-84 zero meridian; units of degrees.
    double m_longitude_deg;
    /// Latitude referenced to WGS-84 equator; units of degrees.
    double m_latitude_deg;
    /// Height referenced to WGS-84 ellipsoid; units of kilometers.
    double height_km;
};

#endif /* GEODETIC_COORD_H */
