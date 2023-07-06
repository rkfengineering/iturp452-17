#ifndef HEIGHT_SCALED_EXCEEDANCE_GRID_H
#define HEIGHT_SCALED_EXCEEDANCE_GRID_H

#include "BinaryFileReader.h"
#include "ExceedanceDataGrid.h"
#include "GeodeticCoord.h"

#include <vector>
#include <memory>

/** Represent a stack of geodetic data grids, each associated with a value in
 * logarithmic scale (i.e. a probability).
 * Implements interpolation method outlined in P.836-6 which requires extra data for 
 * altitude scaling
 * Only for use with total columnar water vapor content and surface water vapor density
 */
class HeightScaledExceedanceDataGrid : public ExceedanceDataGrid{
public:
    HeightScaledExceedanceDataGrid(const BinaryFileReader &source, 
                const std::string &namePattern, const double &resolution,
                const std::vector<ExceedanceLayer> &exceedanceLayerList, 
                const std::shared_ptr<ExceedanceDataGrid> &waterVaporScaleHeightGrid, 
                const std::shared_ptr<DataGrid> &topoHeightGrid)
        : ExceedanceDataGrid(source, namePattern, resolution, exceedanceLayerList), 
                m_waterVaporScaleHeightGrid(waterVaporScaleHeightGrid), 
                m_topoHeightGrid(topoHeightGrid) {}

    /** Interpolate a particular geodetic location on the grid and a depth value.
     *  Scales data to appropriate height using vsch grid and topo grid
     *  Only for use with total columnar water vapor content and surface water vapor density
     * @param location The desired location.
     * @param depth The depth value to interpolate in log-scale.
     * @return The value interpolated from the grid.
     */
    double interpolate2D_withHeight(const GeodeticCoord &location, const double &percentExceedance) const;

private:
    std::shared_ptr<ExceedanceDataGrid> m_waterVaporScaleHeightGrid;
    std::shared_ptr<DataGrid> m_topoHeightGrid;
};

#endif /* HEIGHT_SCALED_EXCEEDANCE_GRID_H */
