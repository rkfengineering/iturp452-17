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

/*
/// @brief 
/// @param hst       Tx antenna height of the smooth-earth surface (amsl) (m)
/// @param hsr       Rx antenna height of the smooth-earth surface (amsl) (m)
/// @param hstd      Tx effective antenna height for the diffraction model (m)
/// @param hsrd      Rx effective antenna height for the diffraction model (m)
/// @param hte       Tx effective antenna height for the ducting/layer reflection model (m)
/// @param hre       Rx effective antenna height for the ducting/layer reflection model (m)
/// @param hm        Terrain roughness parameter (m)
/// @param dlt       Interfering antenna horizon distance (km)
/// @param dlr       Interfered-with antenna horizon distance (km)
/// @param theta_t   Interfering antenna horizon elevation angle (mrad)
/// @param theta_r   Interfered-with antenna horizon elevation angle (mrad)
/// @param theta_tot Angular distance (mrad)
/// @param pathtype  1 = LOS, 2 = transhorizon
struct SmoothEarthResults{
    double hst;
    double hsr;
    double hstd;
    double hsrd;
    double hte;
    double hre;
    double hm;
    double dlt;
    double dlr;
    double theta_t;
    double theta_r;
    double theta_tot;
    int pathtype;
};

/// @brief Smooth-Earth Effective Antenna Heights from Sections 4,5 of Annex 2 of ITU-R P.452.17
/// @param path     Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
/// @param htg      Tx antenna height above ground level (m)
/// @param hrg      Rx antenna height above ground level (m)
/// @param ae       Median effective Earth's radius (km)
/// @param freqGHz  Frequency (GHz)
/// @return         See SmoothEarthResults
SmoothEarthResults smoothEarthHeights(const PathProfile::Path& path, double htg, double hrg, double ae, double freqGHz);
*/
#endif /* PATH_PROFILE_H */