#include "DataLoader.h"
#include "DataStructures.h"

#include <ostream>
#include <sstream>

using namespace Gas;

DataLoader::DataLoader(const BinaryFileReader& binarySrc) 
			: _topographicHeightGridPtr(std::make_shared<DataGrid>(binarySrc, "TOPO_0DOT5", 0.5)),
			_waterVaporScaleHeightGridPtr(std::make_shared<ExceedanceDataGrid>(binarySrc, "VSCH_%1_v4", 1.125, DataStructures::ITU_LAYER_LIST)),
			_columnarCloudWaterVaporContent(binarySrc, "ESAWRED_%1_v4", 1.125, DataStructures::ITU_LAYER_LIST),
			_surfaceWaterVaporDensity(binarySrc, "RHO_%1_v4", 1.125, DataStructures::ITU_LAYER_LIST, 
						_waterVaporScaleHeightGridPtr, 
						_topographicHeightGridPtr),
			_totalWaterVaporContent(binarySrc, "V_%1_v4", 1.125, DataStructures::ITU_LAYER_LIST, 
						_waterVaporScaleHeightGridPtr, 
						_topographicHeightGridPtr) {}

double DataLoader::fetchTopographicHeight_km(const GeodeticCoord& location) const {
	const double TOPOGRAPHIC_HEIGHT_KM = _topographicHeightGridPtr->interpCubic(location);
	return TOPOGRAPHIC_HEIGHT_KM;
}

double DataLoader::fetchTotalColumnarContentCloudWater(const GeodeticCoord& location, const double& exceed) const {
	const double PERCENT_EXCEED = exceed * 100.0;
	if (PERCENT_EXCEED < 0.0 || PERCENT_EXCEED > 100.0) {
		std::ostringstream oStrStream; oStrStream << "ERROR: DataLoader::fetchTotalColumnarContentCloudWater(): " 
					<< "Exceedance is outside of valid range of [0%, 100%]: " 
					<< PERCENT_EXCEED << "%!";
		throw std::domain_error(oStrStream.str());
	}
	const double COLUMNAR_WATER_VAPOR = _columnarCloudWaterVaporContent.interpolate2D(location, PERCENT_EXCEED);
	return COLUMNAR_WATER_VAPOR;
}

double DataLoader::fetchSurfaceWaterVaporDensity(const GeodeticCoord& location, const double& exceed) const {
	const double PERCENT_EXCEED = exceed * 100.0;
	if (PERCENT_EXCEED < 0.0 || PERCENT_EXCEED > 100.0) {
		std::ostringstream oStrStream; oStrStream << "ERROR: ItuModels13::fetchSurfaceWaterVaporDensity(): " 
					<< "Exceedance is outside of valid range of [0%, 100%]: " << PERCENT_EXCEED << "%!" << std::endl;
		throw std::domain_error(oStrStream.str());
	}
	if (location.height_km < 0.0) {
		std::ostringstream oStrStream; oStrStream << "ERROR: ItuModels13::fetchSurfaceWaterVaporDensity(): " 
					<< "Height must be above ground (>= 0 km): " << location.height_km << " km!" << std::endl;
		throw std::domain_error(oStrStream.str());
	}
	const double SURFACE_WATER_VAPOR_DENSITY = _surfaceWaterVaporDensity.interpolate2D_withHeight(location, PERCENT_EXCEED);
	return SURFACE_WATER_VAPOR_DENSITY;
}


double DataLoader::fetchTotalWaterVaporContent(const GeodeticCoord& location, const double& exceed) const {
	const double PERCENT_EXCEED = exceed * 100.0;
	if (PERCENT_EXCEED < 0.0 || PERCENT_EXCEED > 100.0) {
		std::ostringstream oStrStream; oStrStream << "ERROR: DataLoader::fetchTotalWaterVaporContent(): " 
					<< "Exceedance is outside of valid range of [0%, 100%]: " 
					<< PERCENT_EXCEED << "%!" << std::endl;
		throw std::domain_error(oStrStream.str());
	}
	if (location.height_km < 0.0) {
		std::ostringstream oStrStream; oStrStream << "ERROR: DataLoader::fetchTotalWaterVaporContent(): " 
					<< "Height must be above ground (>= 0 km): " 
					<< location.height_km << " km!" << std::endl;
		throw std::domain_error(oStrStream.str());
	}
	const double WATER_VAPOR_CONTENT = _totalWaterVaporContent.interpolate2D_withHeight(location, PERCENT_EXCEED);
	return WATER_VAPOR_CONTENT;
}