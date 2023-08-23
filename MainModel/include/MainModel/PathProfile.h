#ifndef PATH_PROFILE_H
#define PATH_PROFILE_H

#include <vector>
#include <string>

namespace PathProfile{
	/// Specifies zone type of a profile point
	/// First value is set to 1 to require that this value is set intentionally (default value 0 has no meaning).
    /// WARNING the validation data test profile column labels contradict themselves. Assume A1 and zone type 1 both mean coastal land
	enum ZoneType {
		CoastalLand = 1,
		Inland = 2,
		Sea = 3,
	};

    /// @brief          A point along the profile path
    /// @param d_km     Distance (km)
    /// @param h_asl_m  Height above sea level (m)
    /// @param zone     Zone type (if specified)
    struct ProfilePoint{
        ProfilePoint();
        ProfilePoint(double distance_km, double height_asl_m);
        ProfilePoint(double distance_km, double height_asl_m, ZoneType zonetype);
        double d_km;
        double h_asl_m;
        ZoneType zone;
    };

    /// @brief An ordered vector of ProfilePoints that constitute a path
    /// This class adds a constructor to create a path object from a csv file
    class Path : public std::vector<PathProfile::ProfilePoint>{
        public:   
        Path();
        Path(std::string csvPath);

        /// @brief Calculate the fraction of the total path that has the sea zone type
        /// @return Fraction of the path over sea (omega in P452-17)
        double calcFracOverSea() const;

        /// @brief Calculate the time percentage for which refractive index lapse-rates exceeding 100 N-units/km 
        /// can be expected in the first 100m of the lower atmosphere
        /// @param centerLatitude_deg The latitude (deg) of the path center point
        /// @return Time percentage beta0 (%)
        double calcTimePercentBeta0(const double& centerLatitude_deg) const;

        /// @brief Calculate the longest contiguous inland distance in the profile path
        /// @return the longest contiguous inland distance (km)
        double calcLongestContiguousInlandDistance_km() const;
    };
}
#endif /* PATH_PROFILE_H */