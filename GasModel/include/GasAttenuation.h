#ifndef GAS_ATTENUATION_H
#define GAS_ATTENUATION_H

#include "BinaryFileReader.h"
#include "DataLoader.h"
#include "Enumerations.h"
#include "GeodeticCoord.h"

namespace GasAttenuation {
    double performCalculations(const Gas::DataLoader& dataLoader, 
            const GeodeticCoord& location, const double& freq_GHz, const double& elevationAngle_deg, 
            const double& exceedance, const Enumerations::Season& season,
            const bool& useAnnexOne, const bool& useStandardAtmosphere);
}

#endif /* GAS_ATTENUATION_H */