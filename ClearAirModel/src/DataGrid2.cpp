#include <ClearAirModel/DataGrid2.h>
#include <Common/MathHelpers.h>

#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>

//New data loader in constructor. all other functions are exact copies
std::vector<std::vector<double>> DataGrid2::readGridData(const std::string& sourceFilePath) const{
	std::vector<std::vector<double>> fileContents{};

	//move to function that reads raw data and populates datagrid
	std::ifstream file;
	file.open(sourceFilePath);
	std::string line,value;
	if(file.is_open()){
		while(std::getline(file,line)){
			std::vector<double> row{};
			//split line into vector of doubles (tokenizer)
			std::stringstream linestream(line);
			while(std::getline(linestream, value, ' ')){
				//throw out initial empty value (row starts with spaces)
				if(!value.empty()){
					row.push_back(std::stod(value));
				}
			}
			//push back vector to _dataGrid
			if(!row.empty()){
				fileContents.push_back(row);
			}
		}
	}
	return fileContents;
}

// NOTE: Assumes default start/end values for latitude/longitude bounds
DataGrid2::DataGrid2(const std::string& sourceFilePath, const double& resolution_deg,
			const double& beginLat_deg, const double& endLat_deg, 
			const double& beginLon_deg, const double& endLon_deg) 
		: _resolution_deg(resolution_deg), 
		_startLat_deg(beginLat_deg), _endLat_deg(endLat_deg), 
		_startLon_deg(beginLon_deg), _endLon_deg(endLon_deg) {
	// One extra for border
	double EXPECTED_NUM_COLUMNS = std::round(abs(_endLon_deg - _startLon_deg) / _resolution_deg) + 1;
	// Edge case: There is no extra column padding for the default ITU data grid //WARNING this grid has 241 columns. ignore
	/*if (_startLon_deg == 0.0 && _endLon_deg == 360.0) {
		EXPECTED_NUM_COLUMNS = std::round(360.0 / _resolution_deg);
	}*/
	
	// One extra for border
	const double EXPECTED_NUM_ROWS = std::round(abs(_endLat_deg - _startLat_deg) / _resolution_deg) + 1;

	try {
		_dataGrid = readGridData(sourceFilePath);
	}
	catch (std::exception& err) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: DataGrid::DataGrid(): Failed reading source from file \"" 
					<< sourceFilePath << "\": " << err.what();
		throw std::runtime_error(oStrStream.str());
	}

	if (_dataGrid.size() != EXPECTED_NUM_ROWS) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: DataGrid::DataGrid(): In file \"" 
					<< sourceFilePath << "\": expected rows " << EXPECTED_NUM_ROWS << ", but read in " << _dataGrid.size();
		throw std::runtime_error(oStrStream.str());
	}
	if (_dataGrid[0].size() != EXPECTED_NUM_COLUMNS) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: DataGrid::DataGrid(): In file \"" 
					<< sourceFilePath << "\": expected columns " << EXPECTED_NUM_COLUMNS << ", but read in " << _dataGrid[0].size();
		throw std::runtime_error(oStrStream.str());
	}
}





//Other functions are copied verbatim from DataGrid class



std::vector<DataGridHelpers::BoundingBoxGridPoint> DataGrid2::getBoundingBoxList(const GeodeticCoord& location) const {
	const auto BOUNDING_BOX_INTEGER_PAIRS = DataGridHelpers::calculateBoundingBoxIntegerPairs(location, _resolution_deg, 
		_startLat_deg, _endLat_deg, _startLon_deg, _endLon_deg);

	const bool isLatitudeAscending = _endLat_deg > _startLat_deg;

	return DataGridHelpers::calculateBoundingBoxGridPointList(location, BOUNDING_BOX_INTEGER_PAIRS, _resolution_deg, _startLat_deg, _startLon_deg, isLatitudeAscending);
}

double DataGrid2::interpolate2D(const GeodeticCoord& location, const std::vector<double>& customWeightList) const {
	const auto BOUNDING_BOX_INTEGER_PAIRS = DataGridHelpers::calculateBoundingBoxIntegerPairs(location, _resolution_deg,
		_startLat_deg, _endLat_deg, _startLon_deg, _endLon_deg);
	
	const auto LON_COL_NEIGHBOR_PAIR = BOUNDING_BOX_INTEGER_PAIRS.first;
	const auto LAT_ROW_NEIGHBOR_PAIR = BOUNDING_BOX_INTEGER_PAIRS.second;

	// Values at each grid point
	const uint64_t row0Ind = static_cast<uint64_t>(LAT_ROW_NEIGHBOR_PAIR.lowPoint);
	const uint64_t row1Ind = static_cast<uint64_t>(LAT_ROW_NEIGHBOR_PAIR.highPoint);
	const uint64_t col0Ind = static_cast<uint64_t>(LON_COL_NEIGHBOR_PAIR.lowPoint);
	const uint64_t col1Ind = static_cast<uint64_t>(LON_COL_NEIGHBOR_PAIR.highPoint);
	const double value00 = _dataGrid[row0Ind][col0Ind];
	const double value01 = _dataGrid[row0Ind][col1Ind];
	const double value10 = _dataGrid[row1Ind][col0Ind];
	const double value11 = _dataGrid[row1Ind][col1Ind];

	const double INTERP_RESULT = DataGridHelpers::interpolate2D({ value00, value01, value10, value11 }, customWeightList, LAT_ROW_NEIGHBOR_PAIR.weightFactor, LON_COL_NEIGHBOR_PAIR.weightFactor);

	return INTERP_RESULT;
}

double DataGrid2::interpolate2D(const GeodeticCoord& location) const {
	const auto BOUNDING_BOX_INTEGER_PAIRS = DataGridHelpers::calculateBoundingBoxIntegerPairs(location, _resolution_deg,
		_startLat_deg, _endLat_deg, _startLon_deg, _endLon_deg);

	const auto LON_COL_NEIGHBOR_PAIR = BOUNDING_BOX_INTEGER_PAIRS.first;
	const auto LAT_ROW_NEIGHBOR_PAIR = BOUNDING_BOX_INTEGER_PAIRS.second;

	// Values at each grid point
	const uint64_t row0Ind = static_cast<uint64_t>(LAT_ROW_NEIGHBOR_PAIR.lowPoint);
	const uint64_t row1Ind = static_cast<uint64_t>(LAT_ROW_NEIGHBOR_PAIR.highPoint);
	const uint64_t col0Ind = static_cast<uint64_t>(LON_COL_NEIGHBOR_PAIR.lowPoint);
	const uint64_t col1Ind = static_cast<uint64_t>(LON_COL_NEIGHBOR_PAIR.highPoint);
	const double value00 = _dataGrid[row0Ind][col0Ind];
	const double value01 = _dataGrid[row0Ind][col1Ind];
	const double value10 = _dataGrid[row1Ind][col0Ind];
	const double value11 = _dataGrid[row1Ind][col1Ind];

	const double INTERP_RESULT = DataGridHelpers::interpolate2D({ value00, value01, value10, value11 }, LAT_ROW_NEIGHBOR_PAIR.weightFactor, LON_COL_NEIGHBOR_PAIR.weightFactor);

	return INTERP_RESULT;
}

double DataGrid2::interpCubic(const GeodeticCoord& location) const {
	const auto BOUNDING_BOX_INTEGER_PAIRS = DataGridHelpers::calculateBoundingBoxIntegerPairs(location, _resolution_deg,
		_startLat_deg, _endLat_deg, _startLon_deg, _endLon_deg);

	const auto LON_COL_NEIGHBOR_PAIR = BOUNDING_BOX_INTEGER_PAIRS.first;
	const auto LAT_ROW_NEIGHBOR_PAIR = BOUNDING_BOX_INTEGER_PAIRS.second;

	const uint16_t GRID_AXIS_LENGTH = 4;

	std::vector<double> rowIndList;
	rowIndList.reserve(GRID_AXIS_LENGTH);
	std::vector<double> colIndList;
	colIndList.reserve(GRID_AXIS_LENGTH);

	// Now wrap to the file grid (implicitly casting)
	const double NUM_ROWS = static_cast<double>(_dataGrid.size());
	const double NUM_COLS = static_cast<double>(_dataGrid[0].size());
	rowIndList.push_back(MathHelpers::clampValueWithinAxis(LAT_ROW_NEIGHBOR_PAIR.lowPoint - 1.0, NUM_ROWS - 1.0));
	rowIndList.push_back(MathHelpers::clampValueWithinAxis(LAT_ROW_NEIGHBOR_PAIR.lowPoint, NUM_ROWS - 1.0));
	rowIndList.push_back(MathHelpers::clampValueWithinAxis(LAT_ROW_NEIGHBOR_PAIR.highPoint, NUM_ROWS - 1.0));
	rowIndList.push_back(MathHelpers::clampValueWithinAxis(LAT_ROW_NEIGHBOR_PAIR.highPoint + 1.0, NUM_ROWS - 1.0));

	colIndList.push_back(MathHelpers::unwrapValueAroundAxis(LON_COL_NEIGHBOR_PAIR.lowPoint - 1.0, 0.0, NUM_COLS));
	colIndList.push_back(MathHelpers::unwrapValueAroundAxis(LON_COL_NEIGHBOR_PAIR.lowPoint, 0.0, NUM_COLS));
	colIndList.push_back(MathHelpers::unwrapValueAroundAxis(LON_COL_NEIGHBOR_PAIR.highPoint, 0.0, NUM_COLS));
	colIndList.push_back(MathHelpers::unwrapValueAroundAxis(LON_COL_NEIGHBOR_PAIR.highPoint + 1.0, 0.0, NUM_COLS));

	// Extract data values to prepare for interpolation
	std::vector<std::vector<double>> gridValueMatrix;
	gridValueMatrix.reserve(GRID_AXIS_LENGTH);
	uint64_t currentRowInd, currentColInd;
	for (uint16_t rowInd = 0; rowInd < GRID_AXIS_LENGTH; rowInd++) {
		std::vector<double> currentRow;
		currentRow.reserve(GRID_AXIS_LENGTH);
		for (uint16_t colInd = 0; colInd < GRID_AXIS_LENGTH; colInd++) {
			currentRowInd = static_cast<uint64_t>(rowIndList[rowInd]);
			currentColInd = static_cast<uint64_t>(colIndList[colInd]);
			currentRow.push_back(_dataGrid[currentRowInd][currentColInd]);
		}
		gridValueMatrix.push_back(currentRow);
	}
	const double INTERP_RESULT = MathHelpers::calculateBicubicInterpolation(gridValueMatrix, LAT_ROW_NEIGHBOR_PAIR.weightFactor, LON_COL_NEIGHBOR_PAIR.weightFactor);

	return INTERP_RESULT;
}
