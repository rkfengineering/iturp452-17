#include <ITUR_P452/DataLoader.h>

#include <ostream>
#include <sstream>
#include <filesystem>

using namespace ITUR_P452;

//data grid object to store surface refractivity (N0) data
const DataGridTxt DataLoader::_seaLevelSurfaceRefractivityMap{CMAKE_CLEARAIR_SRC_DIR / std::filesystem::path("data/N050.TXT"),1.5};
const DataGridTxt DataLoader::_radioRefractivityIndexLapseRateMap{CMAKE_CLEARAIR_SRC_DIR / std::filesystem::path("data/DN50.TXT"),1.5};

double DataLoader::fetchSeaLevelSurfaceRefractivity(const GeodeticCoord& location){
    return _seaLevelSurfaceRefractivityMap.interpolate2D(location);
}
double DataLoader::fetchRadioRefractivityIndexLapseRate(const GeodeticCoord& location){
    return _radioRefractivityIndexLapseRateMap.interpolate2D(location);
}