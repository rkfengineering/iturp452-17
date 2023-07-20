## Notes

Clutterloss has been written but not tested yet

Diffraction loss has been partially tested
* delta bullington loss (and helper functions) passes.
* diffraction loss for time percentage is not tested yet; the first set of validation data only covers the delta bullington function 

Effective Earth has been partially tested
* AMSL, diffractionModelHeights, horizon distance and elevation passes. Ducting Model heights has not been written yet

The constructor for PathProfile from csv has been tested. The testing is not rigorous, nor is the function particularly well designed
* class function for getting fraction of path over sea has been tested

InvCumNorm math function has been tested. This may be moved to a common library in the future. 

Basic Transmission loss has been tested
* including gas attenuation and multipath focusing corrections
* WARNING the official validation data uses the gas models with frequencies below 1 GHz

TODO - check naming to clarify between meters above sea level and meters above mean sea level. some variables might be mixed up

## Introduction
This library implements most of the algorithms found in ITU-P. 452-17. A few of the algorithms have not been implemented yet and are currently being worked on (as of July 14th, 2023). 
* Implementation is based on: [ITU P.452-17 Recommendation](https://www.itu.int/rec/R-REC-P.452-17-202109-I/en)
* Validation of our implementation is based on: [ITU Official Validation examples for the delta Bullington diffraction prediction method](https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/Validation%20examples%20for%20the%20delta%20Bullington%20diffraction%20prediction%20method.docx)
and also on: [ITU Official Validation examples for software implementation of Recommendation ITU-R P.452](https://www.itu.int/en/ITU-R/study-groups/rsg3/ionotropospheric/R19-WP3M-C-0364!N18-P2!ZIP-E.zip)


