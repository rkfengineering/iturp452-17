#ifndef ITU_MULTI_GRID_HELPERS_H
#define ITU_MULTI_GRID_HELPERS_H

#include "BinaryFileReader.h"
#include "DataGrid.h"

#include <string>
#include <vector>

namespace DataGridLayerHelpers {
	/// @brief Given a raw layer name:
    /// 1) validate that it contains a replacement string (which is "%1")
	/// 2) replace that replacement string with a desired layer ID
	/// @param rawLayerName Raw layer name containing replacement string
	/// @param layerID Layer ID which will replace the replacement string
	/// @return Final validated layer name containing layer ID
	std::string validateAndGenerateLayerName(const std::string &rawLayerName, const std::string &layerID);

    std::vector<DataGrid> generateMonthlyDataGridList(const BinaryFileReader& source, 
				const std::string& rawLayerName, const double& resolution_deg, 
				const double& beginLat_deg, const double& endLat_deg, 
				const double& beginLon_deg, const double& endLon_deg);
} // end namespace DataGridLayerHelpers

#endif /* ITU_MULTI_GRID_HELPERS_H */