## Introduction
This library implements the algorithms found in ITU-P. 452-17 associated with Clear Air predictions. The algorithms for calculating transmission loss have been implemented (as of August 2nd, 2023). 
* Implementation is based on: [ITU P.452-17 Recommendation](https://www.itu.int/rec/R-REC-P.452-17-202109-I/en)
* Validation of our implementation is based on: [ITU Official Validation examples for the delta Bullington diffraction prediction method](https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/Validation%20examples%20for%20the%20delta%20Bullington%20diffraction%20prediction%20method.docx)
and also on: [ITU Official Validation examples for software implementation of Recommendation ITU-R P.452](https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/R19-WP3M-C-0364!N18-P2!ZIP-E.zip)

## How to Use
The loss is calculated by first creating an instance of a p452_TotalAttenuation object. The inputs for the constructor are shown below. 
```
    /// @brief Calculates Basic transmission loss (dB), not exceeded for the required annual percentage time, p,
    /// @param freq_GHz             Frequency (GHz)
    /// @param p_percent            Percentage of time not exceeded (%), 0<p<=50
    /// @param path                 Contains vector of terrain profile distances from Tx (km) and heights (amsl) (m)
    /// @param height_tx_m          Tx Antenna height (m)
    /// @param height_rx_m          Rx Antenna height (m)
    /// @param centerLatitude_deg   The latitude of the path center point (deg) 
    /// @param txHorizonGain_dBi    Tx Antenna directional gain towards the horizon along the path (dB)
    /// @param rxHorizonGain_dBi    Rx Antenna directional gain towards the horizon along the path (dB)
    /// @param pol                  Polarization type (horizontal or vertical)
    /// @param dist_coast_tx_km     Distance over land from Tx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param dist_coast_rx_km     Distance over land from Rx to the coast along the profile path (km) (0 for terminal at sea)
    /// @param deltaN               Average radio-refractive index lapse-rate through the lowest 1km of the atmosphere (positive value) 
    /// @param surfaceRefractivity  Sea Level Surface Refractivity (N0) (N-Units)
    /// @param temp_K               Temperature (K)
    /// @param dryPressure_hPa      Dry air pressure (hPa)
    /// @param tx_clutterType       Clutter Category Type at Tx 
    /// @param rx_clutterType       Clutter Category Type at Rx 
    p452_TotalAttenuation(const double& freq_GHz, const double& p_percent, const PathProfile::Path path, 
            const double& height_tx_m, const double& height_rx_m, const double& centerLatitude_deg, 
            const double& txHorizonGain_dBi, const double& rxHorizonGain_dBi, 
            const Enumerations::PolarizationType& pol, const double& dist_coast_tx_km, 
            const double& dist_coast_rx_km, const double& deltaN, const double& surfaceRefractivity,
            const double& temp_K, const double& dryPressure_hPa, const ClutterType& tx_clutterType, 
            const ClutterType& rx_clutterType);
```
The class function `getTotalTransmissionLoss_dB()` is then used to get the loss value. An example is shown below.
```
double LOSS_VAL = myP452Model.getTotalTransmissionLoss_dB();
```

Path objects can be created from csv files. See the tests/test_paths folder for csv examples.
```    
const PathProfile::Path my_path("my_full_filepath.csv");
```

The following ClutterType values are available under the ClearAirModel namespace:
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

Basic Propagation
* Free Space Propagation Tested
* Free Space Propagation with Gas Attenuation Tested
* Free Space Propagation with Gas Attenuation and Multipath Focusing Correction Tested

Diffraction Loss
* Diffraction loss for p percent of time Tested
* Delta Bullington Loss Tested
* Bullington Loss Tested with smooth path and actual terrain
* Spherical Earth Diffraction Loss Tested
* Spherical Earth Diffraction Loss First Term Tested 
* Spherical Earth Diffraction Loss First Term helper function not directly tested (used for First Term calculations)
* Diffraction Model effective heights Tested

Anomalous Propagation
* Anomalous Propagation Loss Tested
* Ducting Model effective heights Tested
* Terrain Roughness calculation Tested
* Fixed Coupling Loss and Time percentage/angular-distance loss helper functions not directly Tested (used for total Anomalous Propagation Loss)

Tropospheric Scatter
* Tropospheric Scatter Loss Tested
* Fetching data from N050 data map Tested

Clutter/height gain Model
* Height Gain Model antenna heights Tested
* Clutter losses not directly Tested (used in total loss calculation)
* modified path not directly Tested (used in horizon angles and distances calculations)
* No Clutter, Dense Urban, Industrial Zone, Dense Suburban Clutter Types Tested
* Accessing Clutter lookup table not directly tested (used to test different clutter types)

Other Helper Functions
* Inverse Cumulative Normal function Tested
* Frequency to Wavelength conversion Tested
* Median Effective Earth Radius calculation Tested
* Function for finding least squares straight path endpoints Tested (helper function for calculating diffraction model and ducting model effective heights)
* Horizon Angles and Distance calculations Tested
* Gas Attenuation not directly tested (used in Free Space Propagation with Gas Attenuation)
* Path angular distance calculation Tested (Used in Anomolous Propagation and Tropospheric Scatter Model)
* Loading Data Grids from .TXT files Tested
* Fetching data from DN50 data map Tested

Path Profile Functions
* Constructor from csv file Tested (with zone data and no zone data)
* Calculating Fraction of path over sea Tested (only one path with sea data was available)
* Time percentage for which refractive index lapse-rates exceeding 100 N-units/km can be expected in the first 100m of the lower atmosphere Tested with two different center latitudes
* Longest Contiguous Inland Distance Tested

Total Attenuation Class
* Total Attenuation is Tested
* Calculating Slope and Path Blending Interpolation Parameters is not directly tested (used to calculate Total Attenuation)

## Notes

* Clutterloss has been tested with Dense Urban, Industrial Zone, Dense Suburban, and No Clutter conditions
* Total Transmission Loss for the clear air model has been tested on all available validation data
* The constructor for PathProfile from csv has been tested. The testing is not rigorous, nor is the function particularly well designed

* InvCumNorm math function has been tested. This may be moved to a common library in the future. 
* DataGrid2 may be integrated into the DataGrid class from the ituModels common library in the future

* WARNING the official validation data uses the gas models with frequencies below 1 GHz, which is not in the domain recommended by ITU-R P.676

* TODO - check naming to clarify between meters above sea level and meters above mean sea level. some variables might be mixed up
