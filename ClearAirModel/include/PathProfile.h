#ifndef PATH_PROFILE_H
#define PATH_PROFILE_H

#include <vector>
#include <string>

namespace PathProfile{
    /// <summary>
	/// Specifies zone type of a profile point
	/// First value is set to 1 to require that this value is set intentionally (default value 0 has no meaning).
	/// </summary>
	enum ZoneType {
		CoastalLand = 1,
		Inland = 2,
		Sea = 3,
	};

    /// @brief          A point along the profile path
    /// @param d_km     Distance (km)
    /// @param h_masl   Height above sea level (m)
    /// @param zone     Zone type (if specified)
    struct ProfilePoint{
        ProfilePoint();
        ProfilePoint(double distance_km, double height_masl);
        ProfilePoint(double distance_km, double height_masl, ZoneType zonetype);
        double d_km;
        double h_masl;
        ZoneType zone;
    };

    /// @brief An ordered vector of ProfilePoints that constitute a path
    /// This class adds a constructor to create a path object from a csv file
    class Path : public std::vector<PathProfile::ProfilePoint>{
        public:   
        Path();
        Path(std::string csvPath);
    };
}
#endif /* PATH_PROFILE_H */