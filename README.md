## Notes

Clutterloss has been written but not tested yet

Diffraction loss has been partially tested
*bullington loss passes. still need to debug spherical earth, delta bullington, and total diffraction loss functions 

Effective Earth has been partially tested
*AMSL, diffractionModelHeights passes. Ducting Model heights has not been written yet

The constructor for PathProfile from csv has been tested. The testing is not rigorous, nor is the function particularly well designed

InvCumNorm math function has been tested. This may be moved to a common library in the future. 

Template README

## Introduction
This library implements most of the algorithms found in ITU-P. 452-17. A few of the algorithms have not been implemented yet and are currently being worked on (as of July 14th, 2023). 
* Implementation is based on: [ITU P.452-17 Recommendation](https://www.itu.int/rec/R-REC-P.618-13-201712-I/en)
* Validation of our implementation is based on: [ITU Official CG-3M3J-13 Validation Examples](https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/CG-3M3J-13-ValEx-Rev6.1.2.xlsx)

There exist two implementations of this library:
* An older version written in C++ using the MSVC compiler, built on Windows using Visual Studio. This version can be found in `ItuExcel/`.
* A newer version written in C++ using the GCC/G++ compiler, built on Linux (Fedora/Ubuntu). This version can be found in `itu_linux/`.

Additionally, a wrapper has been developed to simplify the input/output structure of this library so that it could be made available as a collection of user-defined Excel functions. That wrapper (and all of its required dependencies) and a sample Excel workbook can be found in `itu_excel_example/`.

## List of all ITU-R Recommendations which are partially or fully implemented in this library

| Identifier  | Official Title |
| ------------- | ------------- |
| ITU-R P.453-13 	| The radio refractive index: its formula and refractivity data  |
| ITU-R P.530-17  	| Propagation data and prediction methods required for the design of terrestrial line-of-sight systems  |
| ITU-R P.618-13 	| Propagation data and prediction methods required for the design of Earth-space telecommunication systems  |
| ITU-R P.676-12  	| Attenuation by atmospheric gases  |
| ITU-R P.835-6 	| Reference Standard Atmospheres  |
| ITU-R P.836-6  	| Water vapour: surface density and total columnar content  |
| ITU-R P.837-7 	| Characteristics of precipitation for propagation modelling  |
| ITU-R P.838-3  	| Specific attenuation model for rain for use in prediction methods |
| ITU-R P.839-4 	| Rain height model for prediction methods  |
| ITU-R P.840-8  	| Attenuation due to clouds and fog |
| ITU-R P.1144-10  	| Interpolation methods for the geophysical properties used to compute propagation effects |
| ITU-R P.1510-1 	| Mean surface temperature |
| ITU-R P.1511-2  	| Topography for Earth-to-space propagation modelling |
| ITU-R P.1623-1  	| Prediction method of fade dynamics on Earth-space paths |
| ITU-R P.1853-1  	| Tropospheric attenuation time series synthesis |

## Functions Defined in Wrapper DLL

* initializeLibrary
	+ `Input:` ituDataDirectoryPath
* interpTopoHeight
	+ `Input:` long_deg, lat_deg
* calculatePathLoss
	+ `Input:` freq_GHz, range_km
* calculateRainHeight
	+ `Input:` long_deg, lat_deg
* calculateRainRate
	+ `Input:` long_deg, lat_deg, exceedance
* calculateRainAttenuation
	+ `Input:` long_deg, lat_deg, height_km, freq_GHz, elDeg, availability, rainHeight_km, rainRate_mmPerHr, polarizationCode
* calculateGasAttenuation_annex1
	+ `Input:` long_deg, lat_deg, height_km, freq_GHz, elevAngle_deg, availability, season
* calculateGasAttenuation_annex2
	+ `Input:` long_deg, lat_deg, height_km, freq_GHz, elevAngle_deg, availability, season
* calculateCrossPolarLoss
    + `Input:` freq_GHz, elevAngle_deg, availability,
	polarizationCode, rainAttenuation_dB
* calculateWetRefractivity
	+ `Input:` long_deg, lat_deg, exceedance
* calculateScintillationAttenuation
	+ `Input:` wetRefractivity, freq_GHz, elevAngle_deg, antDiameter_m, antEfficiency, availability
* calculateCloudAttenuation
	+ `Input:` long_deg, lat_deg, freq_GHz, elevAngle_deg, availability
* calculateTotalAttenuation
	+ `Input:` gasAttenuation_dB, rainAttenuation_dB, cloudAttenuation_dB, scintillationAttenuation_dB

## ITU-R Data from the Recommendations implemented

The library requires a certain set of datafiles to work correctly. The location of these datafiles is required when calling `initializeLibrary` in the wrapper DLL. Those datafiles include:

*  **Topographic heights-**  Data from @xref{ITU-R-P836-6}
*  **Water vapor scale height for intermediate calculations-** Data from @xref{ITU-R-P836-6}
*  **Water vapor density-** Data from @xref{ITU-R-P836-6}
*  **Rain heights-** Data from @xref{ITU-R-P839}.
*  **Monthly mean total rainfall-** Data from @xref{ITU-R-P837-7}
*  **Annual Rain rate exceeded p=0.01%-** Data from @xref{ITU-R-P837-7}
*  **Water vapor density for scintillation-** Data from @xref{ITU-R-P836}
*  **Cloud vapor content-** Data from @xref{ITU-R-P840}
*  **Wet term of the surface refractivity-** Data from @xref{ITU-R-P453-14}
*  **Annual mean surface temperature-** Data from @xref{ITU-R-P1510-1}
*  **Monthly mean surface temperature-** Data from @xref{ITU-R-P1510-1}


**Installation Options**
---

1. Clone code in Visual Studio with [`https://github.com/rkfengineering/ituModels.git`]
2. Clean & Build code by right clicking on
    + `Clean`
    + `Build`

2. Observe final output showing dllMain.cpp creating ItuExcel.dll under
	+ `\itu_excel\x64\Debug\ItuExcel.dll.`  
3. Open `quick rain macro.xlsm` 
    + `Go to developer Tab on Excel`	
    + `Click on Macro`
    + `Edit path so it points to ItuExcel.dll or place folder in the set path in the module`
    + `Edit path for itu_init so that it points to ItuData folder`  	
4. Save and Close `quick rain macro.xlsm` and reopen
5. Wait for start up Msg box to say `ITU rain and gas data loaded` and start computing cells


**Configuration Option Troubleshooting**
---

1. File path not found

    + Check if you have the correct ItuData folder?
        - Sometimes it is possible to set path to ItuData folder with missing binaryGrid files.
    + Check that directory path is set correctly for DLL and ItuData folder?
        - In correct path set will crash Excel or bring you into Debug mode.
    + Check for addtional missing windows dlls?
        - Might me missing necessary windows dlls these can be found in Libraries\Debug
        - Include these missing windows Dlls where ou stored you ItuExcel.dll files
