#ifndef LATLON_COORD_H
#define LATLON_COORD_H

#include <iostream>

class LatLonCoord {
public:
    // Constructor
    LatLonCoord(const double& lat_deg, const double& lon_deg, const double& height_km = 0.0) : 
                lat_deg_(lat_deg), lon_deg_(lon_deg), height_meters_(height_km * 1.0e3) {}

    LatLonCoord() : lat_deg_(0.0), lon_deg_(0.0), height_meters_(0.0) {}

    // Methods
    // operator==
    bool operator== (const LatLonCoord& otherLatLon)
    {
        return lat_deg_ == otherLatLon.lat_deg_ 
                    && lon_deg_ == otherLatLon.lon_deg_
                    && height_meters_ == otherLatLon.height_meters_;
    }

    bool operator== (const LatLonCoord& otherLatLon) const
    {
        return lat_deg_ == otherLatLon.lat_deg_ 
                    && lon_deg_ == otherLatLon.lon_deg_
                    && height_meters_ == otherLatLon.height_meters_;
    }

    // operator!=
    bool operator!= (const LatLonCoord& otherLatLon)
    {
        return lat_deg_ != otherLatLon.lat_deg_ 
                    || lon_deg_ != otherLatLon.lon_deg_
                    || height_meters_ != otherLatLon.height_meters_;
    }

    bool operator!= (const LatLonCoord& otherLatLon) const
    {
        return lat_deg_ != otherLatLon.lat_deg_ 
                    || lon_deg_ != otherLatLon.lon_deg_ 
                    || height_meters_ != otherLatLon.height_meters_;
    }

    // operator<
    bool operator< (const LatLonCoord& otherLatLon)
    {
        if (lat_deg_ == otherLatLon.lat_deg_) {
            if (lon_deg_ == otherLatLon.lon_deg_) {
                return height_meters_ < otherLatLon.height_meters_;
            }
            return lon_deg_ < otherLatLon.lon_deg_;
        }

        return lat_deg_ < otherLatLon.lat_deg_;
    }

    bool operator< (const LatLonCoord& otherLatLon) const
    {
        if (lat_deg_ == otherLatLon.lat_deg_) {
            if (lon_deg_ == otherLatLon.lon_deg_) {
                return height_meters_ < otherLatLon.height_meters_;
            }
            return lon_deg_ < otherLatLon.lon_deg_;
        }

        return lat_deg_ < otherLatLon.lat_deg_;
    }

    // operator>
    bool operator> (const LatLonCoord& otherLatLon)
    {
        if (lat_deg_ == otherLatLon.lat_deg_) {
            if (lon_deg_ == otherLatLon.lon_deg_) {
                return height_meters_ > otherLatLon.height_meters_;
            }
            return lon_deg_ > otherLatLon.lon_deg_;
        }

        return lat_deg_ > otherLatLon.lat_deg_;
    }
    
    bool operator> (const LatLonCoord& otherLatLon) const
    {
        if (lat_deg_ == otherLatLon.lat_deg_) {
            if (lon_deg_ == otherLatLon.lon_deg_) {
                return height_meters_ > otherLatLon.height_meters_;
            }
            return lon_deg_ > otherLatLon.lon_deg_;
        }

        return lat_deg_ > otherLatLon.lat_deg_;
    }

    // operator<=
    bool operator<= (const LatLonCoord& otherLatLon)
    {
        if (lat_deg_ == otherLatLon.lat_deg_) {
            if (lon_deg_ == otherLatLon.lon_deg_) {
                return height_meters_ <= otherLatLon.height_meters_;
            }
            return lon_deg_ <= otherLatLon.lon_deg_;
        }

        return lat_deg_ <= otherLatLon.lat_deg_;
    }
    
    bool operator<= (const LatLonCoord& otherLatLon) const
    {
        if (lat_deg_ == otherLatLon.lat_deg_) {
            if (lon_deg_ == otherLatLon.lon_deg_) {
                return height_meters_ <= otherLatLon.height_meters_;
            }
            return lon_deg_ <= otherLatLon.lon_deg_;
        }

        return lat_deg_ <= otherLatLon.lat_deg_;
    }

    // operator>=
    bool operator>= (const LatLonCoord& otherLatLon)
    {
        if (lat_deg_ == otherLatLon.lat_deg_) {
            if (lon_deg_ == otherLatLon.lon_deg_) {
                return height_meters_ >= otherLatLon.height_meters_;
            }
            return lon_deg_ >= otherLatLon.lon_deg_;
        }

        return lat_deg_ >= otherLatLon.lat_deg_;
    }

    bool operator>= (const LatLonCoord& otherLatLon) const
    {
        if (lat_deg_ == otherLatLon.lat_deg_) {
            if (lon_deg_ == otherLatLon.lon_deg_) {
                return height_meters_ >= otherLatLon.height_meters_;
            }
            return lon_deg_ >= otherLatLon.lon_deg_;
        }

        return lat_deg_ >= otherLatLon.lat_deg_;
    }

    friend std::ostream& operator<< (std::ostream &out, LatLonCoord const& data) {
        out << "(Lat = " << data.lat_deg_ << " deg, Lon = " << data.lon_deg_ << " deg)";
        return out;
    }

    // Properties
    double lat_deg_;
    double lon_deg_;
    double height_meters_;
};

#endif // LATLON_COORD_H