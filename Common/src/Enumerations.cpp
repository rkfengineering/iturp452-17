#include "Enumerations.h"

const double Enumerations::determineAngleOfPolarization_deg(const PolarizationType& polarType) {
    switch (polarType) {
    case PolarizationType::HorizontalPolarized:
        return 0.0;
    case PolarizationType::CircularPolarized:
        return 45.0;
    case PolarizationType::VerticalPolarized:
        return 90.0;
    default:
        std::ostringstream oStrStream;
        oStrStream << "ERROR: Enumerations::determineAngleOfPolarization_deg(): The polarization type given is invalid: " << int(polarType) << "!";
        throw std::domain_error(oStrStream.str());
    }
}

void Enumerations::validateSeason(const Season& season) {
    switch (season) {
    case SummerTime:
    case WinterTime:
        break;
    default:
        std::ostringstream oStrStream;
        oStrStream << "ERROR: Enumerations::validateSeason(): The season given is invalid: " << int(season) << "!";
        throw std::domain_error(oStrStream.str());
    }
}