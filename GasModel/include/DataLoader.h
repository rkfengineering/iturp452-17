#ifndef GAS_DATA_LOADER_H
#define GAS_DATA_LOADER_H

#include "BinaryFileReader.h"
#include "ExceedanceDataGrid.h"
#include "GeodeticCoord.h"
#include "HeightScaledExceedanceDataGrid.h"

#include <memory>

namespace Gas {
    class DataLoader {
    public:
        /// @brief Populates the datasets required to calculate cloud attenuation
        /// @param binarySrc Helper to read .bingrd files containing the datasets
        DataLoader(const BinaryFileReader& binarySrc);

        /// @brief Fetch the topographic height above mean sea level
        /// @param location Desired location (lat/lon)
        /// This data is fetched from a binary grid file taken from ITU-R P.836-6
        double fetchTopographicHeight_km(const GeodeticCoord& location) const;

        /// @brief Fetch the total columnar content of reduced cloud liquid water for a given (lat/lon) location
        /// This data is fetched from a binary grid file taken from ITU-R P.840-8
        /// @param location Desired location (lat/lon)
        /// @param exceed Desired exceedance of the system (0-1)
        /// @return Total columnar content of reduced cloud liquid water (kg/m^2)
        double fetchTotalColumnarContentCloudWater(const GeodeticCoord& location, const double& exceed) const;

        double fetchSurfaceWaterVaporDensity(const GeodeticCoord& location, const double& exceed) const;

        /// @brief Fetch the total water vapor content for a given (lat/lon) location
        /// This data is fetched from a binary grid file taken from ITU-R P.836-6
        /// @param location Desired location (lat/lon)
        /// @param exceed Desired exceedance (0-1)
        /// @return Total water vapor content (kg/m^2)
        double fetchTotalWaterVaporContent(const GeodeticCoord& location, const double& exceed) const;

    private:
        /// Topographic heights. Data from ITU-R P.836-6
        std::shared_ptr<DataGrid> _topographicHeightGridPtr;
        /// Water vapor scale height for intermediate calculations. Data from ITU-R P.836-6
        std::shared_ptr<ExceedanceDataGrid> _waterVaporScaleHeightGridPtr;

        /// Annual values of total columnar content of reduced cloud liquid water, Lred (kg/m2). Data from ITU-R P.840-8
        ExceedanceDataGrid _columnarCloudWaterVaporContent;
        /// Water vapor density for scintillation. Data from ITU-R P.836-6
        HeightScaledExceedanceDataGrid _surfaceWaterVaporDensity;
        /// Total water vapour content, expressed in kg/m2 or, equivalently, in mm of precipitable water. Data from ITU-R P.836-6
        /// NOTE: The annual values of total columnar water vapour content, V, in kg/m2 and corresponding water vapour scale height, vsch, are integral parts of this data
        HeightScaledExceedanceDataGrid _totalWaterVaporContent;
    };
} // end namespace Gas

#endif /* GAS_DATA_LOADER_H */