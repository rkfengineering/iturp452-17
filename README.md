## Introduction
This library implements the algorithms found in ITU-P. 452-17 associated with Clear Air predictions. The algorithms for calculating transmission loss have been implemented (as of August 2nd, 2023). 
* Implementation is based on: [ITU P.452-17 Recommendation](https://www.itu.int/rec/R-REC-P.452-17-202109-I/en)
* Validation of our implementation is based on: [ITU Official Validation examples for the delta Bullington diffraction prediction method](https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/Validation%20examples%20for%20the%20delta%20Bullington%20diffraction%20prediction%20method.docx)
and also on: [ITU Official Validation examples for software implementation of Recommendation ITU-R P.452](https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/R19-WP3M-C-0364!N18-P2!ZIP-E.zip)

## How to Use (User Interface Functions)
The following functions can be found in P452/P452.h:

To calculate path loss using a vector of elevation values, use:
```
double calculateP452Loss_dB(const double& txHeight_m, const double& rxHeight_m, 
        const std::vector<double>& elevationList_m, const double& stepDistance_km, 
        const double& midpoint_lat_deg, const double& midpoint_lon_deg,
        const double& freq_GHz, const double& timePercent, const int& polariz=0,
        const double& txHorizonGain_dBi=0, const double& rxHorizonGain_dBi=0,
        const ClutterModel::ClutterType& txClutterType=ClutterModel::ClutterType::NoClutter,
        const ClutterModel::ClutterType& rxClutterType=ClutterModel::ClutterType::NoClutter);

```
To calculate path loss using terrain map data already loaded into a gdal raster processor, use:
```
double calculateP452Loss_dB(const double& txHeight_m, const double& rxHeight_m, 
        const LatLonCoord startCoord, const LatLonCoord endCoord,
        const std::vector<GdalRasterProcessor>& rasterProcessorList,
        const double& freq_GHz, const double& timePercent, const int& polariz=0,
        const double& txHorizonGain_dBi=0, const double& rxHorizonGain_dBi=0,
        const ClutterModel::ClutterType& txClutterType=ClutterModel::ClutterType::NoClutter,
        const ClutterModel::ClutterType& rxClutterType=ClutterModel::ClutterType::NoClutter);
```
To calculate path loss using terrain map data already loaded into a gdal raster processor
and zone classification using land border data loaded into a gdal vector processor, use:
```
double calculateP452Loss_dB(const double& txHeight_m, const double& rxHeight_m, 
        const LatLonCoord startCoord, const LatLonCoord endCoord,
        const std::vector<GdalRasterProcessor>& rasterProcessorList,
        const GdalVectorProcessor& landBorders,
        const double& freq_GHz, const double& timePercent, const int& polariz=0,
        const double& txHorizonGain_dBi=0, const double& rxHorizonGain_dBi=0,
        const ClutterModel::ClutterType& txClutterType=ClutterModel::ClutterType::NoClutter,
        const ClutterModel::ClutterType& rxClutterType=ClutterModel::ClutterType::NoClutter);
```

## How to Use (Internal Core Algorithm Function)
The following function can be found in MainModel/P452TotalAttenuation.h:

The loss is calculated by first creating an instance of a TotalClearAirAttenuation object. The inputs for the constructor are shown below. 
The output is the basic transmission loss (dB), not exceeded for the required annual percentage time, p.

|Input Parameter | Parameter Description|
|----------------|----------------------|
|freq_GHz             |Frequency (GHz)|
|p_percent            |Required time percentage for which the calculated loss is not exceeded, 0< p< =50|
|path                 |Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)|
|height_tx_m          |Tx Antenna height (m)|
|height_rx_m          |Rx Antenna height (m)|
|centerLatitude_deg   |The latitude of the path center point (deg) |
|txHorizonGain_dBi    |Tx Antenna directional gain towards the horizon along the path (dB)|
|rxHorizonGain_dBi    |Rx Antenna directional gain towards the horizon along the path (dB)|
|pol                  |Polarization type (horizontal or vertical)|
|dist_coast_tx_km     |Distance over land from Tx to the coast along the profile path (km) (0 for terminal at sea)|
|dist_coast_rx_km     |Distance over land from Rx to the coast along the profile path (km) (0 for terminal at sea)|
|deltaN               |Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (positive value) |
|surfaceRefractivity  |Sea Level Surface Refractivity (N0) (N-Units)|
|temp_K               |Temperature (K)|
|dryPressure_hPa      |Dry air pressure (hPa)|
|tx_clutterType       |Clutter Category Type at Tx |
|rx_clutterType       |Clutter Category Type at Rx |

The class function `calcTotalClearAirAttenuation` is then used to get the loss value. An example is shown below.
```
double LOSS_VAL = myP452Model.calcTotalClearAirAttenuation;
```

Path objects can be created from csv files. See the tests/test_paths folder for csv examples.
```    
const PathProfile::Path my_path("my_full_filepath.csv");
```

The following ClutterType values are available under the ClutterModel namespace:
```
enum ClutterType {
    NoClutter = 0,
    HighCropFields,
    ParkLand,
    IrregularlySpacedSparseTrees,
    Orchard_RegularlySpaced,
    SparseHouses,
    VillageCentre,
    DeciduousTrees_IrregularlySpaced,
    DeciduousTrees_RegularlySpaced,
    MixedTreeForest,
    ConiferousTrees_IrregularlySpaced,
    ConiferousTrees_RegularlySpaced,
    TropicalRainForest,
    Suburban,
    DenseSuburban,
    Urban,
    DenseUrban,
    HighRiseUrban,
    IndustrialZone
};
```

## Submodel Tests

|Basic Propagation||
|---|---|
|Free Space Propagation |Tested|
|Free Space Propagation with Gas Attenuation |Tested|
|Free Space Propagation with Gas Attenuation and Multipath Focusing Correction |Tested|

|Diffraction Loss||
|---|---|
|Diffraction loss for p percent of time |Tested|
|Delta Bullington Loss |Tested|
|Bullington Loss |Tested with smooth path and actual terrain|
|Spherical Earth Diffraction Loss |Tested|
|Spherical Earth Diffraction Loss First Term |Tested| 
|Spherical Earth Diffraction Loss First Term single zone function |Not directly tested (used for First Term calculations)|
|Diffraction Model effective heights |Tested|

|Anomalous Propagation||
|---|---|
|Anomalous Propagation Loss |Tested|
|Ducting Model effective heights |Tested|
|Terrain Roughness calculation |Tested|
|Fixed Coupling Loss and Time percentage/angular-distance loss helper functions |Not directly Tested (used for total Anomalous Propagation Loss)|
|Gaseous Attenuation |Not directly Tested (used for total Anomalous Propagation Loss)|

|Tropospheric Scatter||
|---|---|
|Tropospheric Scatter Loss |Tested|
|Fetching data from N050 data map |Tested|

|Clutter/height gain Model||
|---|---|
|Height Gain Model antenna heights |Tested|
|Clutter losses |Not directly Tested (used in total loss calculation)|
|modified path |Not directly Tested (used in horizon angles and distances calculations)|
|Accessing Clutter lookup table| Not directly tested (used to test No Clutter, Dense Urban, Industrial Zone, Dense Suburban)|

|Other Helper Functions||
|---|---|
|Inverse Cumulative Normal function |Tested|
|Frequency to Wavelength conversion |Tested|
|Median Effective Earth Radius calculation |Tested|
|Function for finding least squares straight path endpoints (helper function for calculating diffraction model and ducting model effective heights) |Tested |
|Horizon Angles and Distance calculations |Tested|
|Gas Attenuation |Not directly tested (used in Free Space Propagation with Gas Attenuation)| 
|Path angular distance calculation (Used in Anomolous Propagation and Tropospheric Scatter Model)|Tested |
|Loading Data Grids from .TXT files |Tested|
|Fetching data from DN50 data map |Tested|

|Path Profile Functions||
|---|---|
|Constructor from csv file |Tested (with zone data and no zone data)|
|Calculating Fraction of path over sea |Tested (only one path with sea data was available)|
|Time percentage for which refractive index lapse-rates exceeding 100 N-units/km can be expected in the first 100m of the lower atmosphere |Tested with two different center latitudes|
|Longest Contiguous Inland Distance |Tested|

|Total Clear Air Attenuation Class||
|---|---|
|Total Clear Air Attenuation |Tested|
|Calculating Slope and Path Blending Interpolation Parameters is |Not directly tested (used to calculate Total Clear Air Attenuation)|

## Notes

* Clutterloss has been tested with Dense Urban, Industrial Zone, Dense Suburban, and No Clutter conditions
* Total Transmission Loss for the clear air model has been tested on all available validation data
* The constructor for PathProfile from csv has been tested. The testing is not rigorous, nor is the function particularly well designed

* InvCumNorm math function has been tested. This may be moved to a common library in the future. 
* DataGridTxt may be integrated into the DataGrid class from the ituModels common library in the future

* WARNING the official validation data uses the gas models with frequencies below 1 GHz, which is not in the domain recommended by ITU-R P.676

* TODO - check naming to clarify between meters above sea level and meters above mean sea level. some variables might be mixed up
