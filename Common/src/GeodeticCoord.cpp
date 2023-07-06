#include "GeodeticCoord.h"

#include <cmath>

bool GeodeticCoord::isNull() const{
    return (
        std::isnan(m_longitude_deg) ||
        std::isnan(m_latitude_deg) ||
        std::isnan(height_km)
    );
}