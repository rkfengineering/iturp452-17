#ifndef PROFILE_PATH_H
#define PROFILE_PATH_H

#include <vector>
#include <string>

/// @brief 
/// @param d        Vector of distances di of the i-th profile point (km)
/// @param h        Vector of heights hi of the i-th profile point (masl)
/// @param zone     Vector of Zone types of the i-th profile point: Coastal Land (1), Inland (2), Sea (3)
/// @param length   Number of profile points
class ProfilePath{
    public:   
    ProfilePath();
    ProfilePath(std::string csvPath);
    std::vector<double> d;
    std::vector<double> h;
    std::vector<int> zone;
    int length;
};

#endif /* PROFILE_PATH_H */