#ifndef EXCEEDANCE_DATA_GRID_H
#define EXCEEDANCE_DATA_GRID_H

#include "BinaryFileReader.h"
#include "DataGridHelpers.h"
#include "DataGridLayerHelpers.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

/** Represent a stack of geodetic data grids, each associated with a value in
 * logarithmic scale (i.e. a probability).
 */
class ExceedanceDataGrid{
public:
    /// Identify each grid in the stack
    struct ExceedanceLayer {
        // Create a new layer identity
        ExceedanceLayer(const std::string &layerId, const double &percentExceedance) 
                    : m_layerId(layerId), m_percentExceedance(percentExceedance), 
                    m_logPercentExceedance(std::log(percentExceedance)) {}

        // Identifier used for the file name
        std::string m_layerId;
        // Linear value used to select the layer
        double m_percentExceedance;
        // Logarithmic value used to select the layer
        double m_logPercentExceedance;
    };

    /** Load a set of geodetic grid files all with a common file name format.
     * NOTE: These are assumed to be sorted in ascending order of value.
     * @param source The mapping from names to numeric data.
     * @param commonFileNamePattern The common file name format including the formating
     * string %1, which will be substituted with the three-digit index of
     * each file.
     * @param resolution The grid resolution in degrees (longitude and latitude)
     * @param exceedanceLayerList List of ExceedanceGridLayers which must be populated
     */
    ExceedanceDataGrid(const BinaryFileReader& source, 
            const std::string& commonFileNamePattern, const double& resolution, 
            const std::vector<ExceedanceLayer>& exceedanceLayerList);

    // Free resources
    ~ExceedanceDataGrid() = default;
    
    /** Get the set of grid points nearest to a location.
     *
     * @param location The desired location.
     * @return The nearest grid points weighted by distance.
     */
    std::vector<DataGridHelpers::BoundingBoxGridPoint> getFirstLayerBoundingBox(const GeodeticCoord& location) const;

    /** Interpolate a particular geodetic location on the grid and a depth value.
     * @param location The desired location.
     * @param depth The depth value to interpolate in log-scale.
     * @return The value interpolated from the grid.
     */
    double interpolate2D(const GeodeticCoord &location, const double &percentExceedance) const;
    
protected:
    // DataGrid layers sorted in ascending order by exceedance
    std::map<double, std::shared_ptr<DataGrid>> _logExceedanceToGridMap;
};

#endif /* EXCEEDANCE_DATA_GRID_H */
