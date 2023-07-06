#include "HeightScaledExceedanceDataGrid.h"

#include "MathHelpers.h"

double HeightScaledExceedanceDataGrid::interpolate2D_withHeight(
            const GeodeticCoord &location, const double &percentExceedance) const {
    const double logPercentExceedance = std::log(percentExceedance);
    const std::vector<DataGridHelpers::BoundingBoxGridPoint> firstLayerBoundingBox = getFirstLayerBoundingBox(location);

    std::vector<double> layerWeightList = { 1,1,1,1 };

    //assume nearpoints are in the correct order
    for (uint16_t pointInd = 0; pointInd < firstLayerBoundingBox.size(); pointInd++) {
        const GeodeticCoord currentBoundingBoxCoord = firstLayerBoundingBox[pointInd].m_location;
        const double currentElevation_km = m_topoHeightGrid->interpCubic(currentBoundingBoxCoord);

        layerWeightList[pointInd] = std::exp(-(location.height_km - currentElevation_km)
            / m_waterVaporScaleHeightGrid->interpolate2D(currentBoundingBoxCoord, percentExceedance));
    }

    // If the desired percent exceedance is smaller than the layer with smallest percent exceedance, use that layer (don't extrapolate) 
    if(logPercentExceedance <=  _logExceedanceToGridMap.begin()->first){
        return _logExceedanceToGridMap.begin()->second->interpolate2D(location, layerWeightList);
    }
    // If the desired percent exceedance is larger than the layer with larger percent exceedance, use that layer (don't extrapolate)
    else if(logPercentExceedance >= _logExceedanceToGridMap.rbegin()->first){
        return _logExceedanceToGridMap.rbegin()->second->interpolate2D(location, layerWeightList);
    }
    // Interpolate between the two layers surrounding the desired percent exceedance
    else{
        // First layer at a depth >= current depth
        const auto highLayerIter = _logExceedanceToGridMap.lower_bound(logPercentExceedance);
        // Layer preceeding the higher layer
        const auto lowLayerIter = std::prev(highLayerIter);

        std::vector<double> lowLayerWeightsList = {1,1,1,1};
        std::vector<double> highLayerWeightsList = {1,1,1,1};

        // Calculate weights for each point in the bounding box relative to the desired location for both layers
        for (uint16_t pointInd = 0; pointInd < firstLayerBoundingBox.size(); pointInd++) {
            const GeodeticCoord currentBoundingBoxCoord = firstLayerBoundingBox[pointInd].m_location;
            const double currentElevation_km = m_topoHeightGrid->interpCubic(currentBoundingBoxCoord);
            
            const double ELEVATION_HEIGHT_SCALE = -(location.height_km - currentElevation_km);
            const double LOW_WATER_VAPOR_HEIGHT_SCALE = m_waterVaporScaleHeightGrid->interpolate2D(currentBoundingBoxCoord, std::exp(lowLayerIter->first));
            const double HIGH_WATER_VAPOR_HEIGHT_SCALE = m_waterVaporScaleHeightGrid->interpolate2D(currentBoundingBoxCoord, std::exp(highLayerIter->first));

            lowLayerWeightsList[pointInd] = std::exp(ELEVATION_HEIGHT_SCALE / LOW_WATER_VAPOR_HEIGHT_SCALE);
            highLayerWeightsList[pointInd] = std::exp(ELEVATION_HEIGHT_SCALE / HIGH_WATER_VAPOR_HEIGHT_SCALE);
        }

        // Interpolate based on distance of desired location from the bounding box points
        const double INTERP_LOW_LAYER_VALUE = lowLayerIter->second->interpolate2D(location, lowLayerWeightsList);
        const double INTERP_HIGH_LAYER_VALUE = highLayerIter->second->interpolate2D(location, highLayerWeightsList);

        // Interpolate based on distance of desired exceedance from the surrounding layer's exceedances
        const double HIGH_LAYER_EXCEEDANCE_WEIGHT = (logPercentExceedance - lowLayerIter->first) / (highLayerIter->first - lowLayerIter->first);
        const double FINAL_INTERP_VALUE = MathHelpers::interpolate1D(INTERP_LOW_LAYER_VALUE, INTERP_HIGH_LAYER_VALUE, HIGH_LAYER_EXCEEDANCE_WEIGHT);

        return FINAL_INTERP_VALUE;
    }
}