## Notes

Clutterloss has been tested with Dense Urban, Industrial Zone, Dense Suburban, and No Clutter conditions

Diffraction loss has been fully tested
* calcSphericalEarthDiffraction_firstTerm_helper_dB is indirectly tested through calcSphericalEarthDiffraction_firstTerm_dB
* Diffraction Loss has been tested on fully inland paths as well as mixed terrain paths

The constructor for PathProfile from csv has been tested. The testing is not rigorous, nor is the function particularly well designed
* class function for getting fraction of path over sea has been tested on mixed terrain and fully inland paths
* class function for getting b0 percentage has been tested
* class function for getting longest contiguous inland distance has been tested

InvCumNorm math function has been tested. This may be moved to a common library in the future. 

Basic Transmission loss has been tested
* including gas attenuation and multipath focusing corrections
* WARNING the official validation data uses the gas models with frequencies below 1 GHz

Transmission Loss for anomolous propagation has been tested for LOS and transhorizon paths
Transmission Loss for tropospheric scatter has been tested
Total Transmission Loss for the clear air model has been tested on all available validation data

TODO - check naming to clarify between meters above sea level and meters above mean sea level. some variables might be mixed up

## Introduction
This library implements the algorithms found in ITU-P. 452-17 associated with Clear Air predictions. The algorithms for calculating transmission loss have been implemented (as of August 2nd, 2023). 
* Implementation is based on: [ITU P.452-17 Recommendation](https://www.itu.int/rec/R-REC-P.452-17-202109-I/en)
* Validation of our implementation is based on: [ITU Official Validation examples for the delta Bullington diffraction prediction method](https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/Validation%20examples%20for%20the%20delta%20Bullington%20diffraction%20prediction%20method.docx)
and also on: [ITU Official Validation examples for software implementation of Recommendation ITU-R P.452](https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/R19-WP3M-C-0364!N18-P2!ZIP-E.zip)


