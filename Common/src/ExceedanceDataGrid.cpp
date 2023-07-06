#include "ExceedanceDataGrid.h"

#include "DataGrid.h"
#include "MathHelpers.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

/* ExceedanceDataGrid definitions */
ExceedanceDataGrid::ExceedanceDataGrid(const BinaryFileReader& source, 
            const std::string& commonFileNamePattern, const double& resolution, 
            const std::vector<ExceedanceLayer>& exceedanceLayerList) {
    
    for (auto& exceedanceLayer : exceedanceLayerList) {
        std::string layerNameWithId = DataGridLayerHelpers::validateAndGenerateLayerName(commonFileNamePattern, exceedanceLayer.m_layerId);

        try {
            auto dataGridPtr = std::make_shared<DataGrid>(source, layerNameWithId, resolution);
            _logExceedanceToGridMap.insert(std::make_pair(exceedanceLayer.m_logPercentExceedance, dataGridPtr));
        }
        catch (std::exception& err) {
            std::ostringstream oStrStream;
            oStrStream << "ERROR: ExceedanceDataGrid::ExceedanceDataGrid(): " 
                        << "Something went wrong while processing the current layer (" << layerNameWithId << "): " << err.what();
            throw std::domain_error(oStrStream.str());
        }
    }
}

std::vector<DataGridHelpers::BoundingBoxGridPoint> ExceedanceDataGrid::getFirstLayerBoundingBox(const GeodeticCoord &location) const{
    if(_logExceedanceToGridMap.empty()){
        throw std::logic_error("ERROR: ExceedanceDataGrid::getFirstLayerBoundingBox(): No layers present!");
    }
    return _logExceedanceToGridMap.begin()->second->getBoundingBoxList(location);
}

double ExceedanceDataGrid::interpolate2D(const GeodeticCoord &location, const double& percentExceedance) const{
    const double logPercentExceedance = std::log(percentExceedance);
    
    // If the desired percent exceedance is smaller than the layer with smallest percent exceedance, use that layer (don't extrapolate) 
    if (logPercentExceedance <= _logExceedanceToGridMap.begin()->first) {
        return _logExceedanceToGridMap.begin()->second->interpolate2D(location);
    }
    // If the desired percent exceedance is larger than the layer with larger percent exceedance, use that layer (don't extrapolate)
    else if (logPercentExceedance >= _logExceedanceToGridMap.rbegin()->first) {
        return _logExceedanceToGridMap.rbegin()->second->interpolate2D(location);
    }
    // Interpolate between the two layers surrounding the desired percent exceedance
    else {
        // First layer at a exceedance >= desired exceedance
        const auto highLayerIter = _logExceedanceToGridMap.lower_bound(logPercentExceedance);
        // Layer preceeding the higher layer
        const auto lowLayerIter = std::prev(highLayerIter);

        // Interpolate based on distance of desired exceedance from the surrounding layer's exceedances
        const double HIGH_LAYER_EXCEEDANCE_WEIGHT = (logPercentExceedance - lowLayerIter->first) 
                    / (highLayerIter->first - lowLayerIter->first);
        const double LOW_LAYER_INTERP = lowLayerIter->second->interpolate2D(location);
        const double HIGH_LAYER_INTERP = highLayerIter->second->interpolate2D(location);

        if (HIGH_LAYER_EXCEEDANCE_WEIGHT == 1.0) {
            return HIGH_LAYER_INTERP;
        }
        else if (HIGH_LAYER_EXCEEDANCE_WEIGHT == 0.0) {
            return LOW_LAYER_INTERP;
        }
        else {
            return MathHelpers::interpolate1D(LOW_LAYER_INTERP, HIGH_LAYER_INTERP, HIGH_LAYER_EXCEEDANCE_WEIGHT);
        }
    }
}
