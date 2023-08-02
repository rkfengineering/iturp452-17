#ifndef DATA_GRID_2_H
#define DATA_GRID_2_H

#include <Common/DataGridHelpers.h>
#include <Common/GeodeticCoord.h>

#include <cstdint>
#include <string>
#include <vector>

/// <summary>
/// TODO integrate with DataGrid class. Takes txt source files instead of .bingrd files
// WARNING not in a namespace
/// </summary>
class DataGrid2 {
public:

    /// @brief Load a grid of data from a binary file which provides world-wide coverage for that data
    /// NOTE: Most ITU data grids are defined from 0->360 deg longitude (not including 360 deg) and -90->90 deg latitude (including 90 deg)
    /// @param sourceFilePath Filepath to the .txt source data
    /// @param resolution_deg Resolution of steps between rows or columns in the data grid (deg)
    /// @param beginLat_deg Starting latitude of the data grid (defaults to 90) (deg)
    /// @param endLat_deg Ending latitude of the data grid (defaults to -90) (deg)
    /// @param beginLon_deg Starting longitude of the data grid (defaults to 0) (deg)
    /// @param endLon_deg Ending longitude of the dat agrid (defaults to 360) (deg)
    DataGrid2(const std::string& sourceFilePath, const double& resolution_deg,
                const double& beginLat_deg = 90.0, const double& endLat_deg = -90.0, 
                const double& beginLon_deg = 0.0, const double& endLon_deg = 360.0);

    std::vector<std::vector<double>> readGridData(const std::string& sourceFilePath) const;

    /// @brief Calculates the bounding-box of coordinates on this data grid which contain the given location
    /// @param location Location which must be bounded within this data grid
    /// @return Bounding-box containing the given location
    std::vector<DataGridHelpers::BoundingBoxGridPoint> getBoundingBoxList(const GeodeticCoord& location) const;

    /// @brief Bi-linear interpolation for a value at the given location from a 2D _data matrix (4-point interpolation)
    /// NOTE: Result is weighted by applying the given weights to the bounding box values
    /// @param location Location where a value is needed
    /// @param customWeightList Weights to be applied to the values at coordinates (0,0), (0,1), (1,0), and (1,1), respectively
    /// @return Custom-weighted interpolation result
    double interpolate2D(const GeodeticCoord& location, const std::vector<double>& customWeightList) const;

    /// @brief Bi-linear interpolation for a value at the given location from a 2D _data matrix (4-point interpolation)
    /// NOTE: Result is weighted by applying the given weights to the bounding box values
    /// @param location Location where a value is needed
    /// @return Interpolation result at given location
    double interpolate2D(const GeodeticCoord& location) const;

    /// @brief Bi-cubic interpolation for a value at the given location from a 2D _data matrix (16-point interpolation)
    /// @param location Location where a value is needed
    /// @return Interpolation result at given location
    double interpCubic(const GeodeticCoord& location) const;

private:
    /// Spacing between points in _dataGrid
    double _resolution_deg;

    /// Boundaries of the data grid in lat/lon space
    double _startLat_deg;
    double _endLat_deg;
    double _startLon_deg;
    double _endLon_deg;

    /// 2D matrix containing data
    /// Rows = latitude, Columns = longitude
    std::vector<std::vector<double>> _dataGrid;
};

#endif /* DATA_GRID_2_H */