#include "DataGridLayerHelpers.h"
#include "DataGrid.h"

#include <iomanip>

std::string DataGridLayerHelpers::validateAndGenerateLayerName(const std::string &rawLayerName, const std::string &layerID)
{
    // Find the expected "%1" symbol in the layer name
    const size_t idPosition = rawLayerName.find("%1");
    if (idPosition == std::string::npos) {
        std::ostringstream oStrStream;
        oStrStream << "ERROR: DataGridLayerHelpers::validateAndGenerateLayerName(): The given layer name did not contain the expected \"%1\" symbol: " << rawLayerName << "!";
        throw std::domain_error(oStrStream.str());
    }

    // Replace the "%1" symbol in the layer name with the desired layer ID
    std::string layerNameWithId(rawLayerName);
    layerNameWithId.replace(idPosition, 2, layerID);

    return layerNameWithId;
}

std::vector<DataGrid> DataGridLayerHelpers::generateMonthlyDataGridList(const BinaryFileReader& source, 
            const std::string& rawLayerName, const double& resolution_deg, 
            const double& beginLat_deg, const double& endLat_deg, 
            const double& beginLon_deg, const double& endLon_deg) {
    std::vector<DataGrid> monthlyDataGridList;
    for (uint16_t monthInd = 1; monthInd <= 12; monthInd++) {
        // Generate month-based layer ID
        std::ostringstream monthIdStream;
        monthIdStream << std::setfill('0') << std::setw(2) << monthInd;
        std::string layerNameWithMonth = validateAndGenerateLayerName(rawLayerName, monthIdStream.str());

        // Populate layer list with monthly data grids in order from first to last month index
        try {
            const DataGrid monthlyDataGrid(source, layerNameWithMonth, resolution_deg, beginLat_deg, endLat_deg, beginLon_deg, endLon_deg);
            monthlyDataGridList.push_back(monthlyDataGrid);
        }
        catch (std::exception& err) {
            std::ostringstream oStrStream;
            oStrStream << "ERROR: DataGridLayerHelpers::generateMonthlyDataGridList(): " 
                        << "Something went wrong will processing the current month's layer (" 
                        << layerNameWithMonth << "): " << err.what();
            throw std::domain_error(oStrStream.str());
        }
    }

    return monthlyDataGridList;
}