#include "GasAttenuationHelpers.h"

namespace {
    double constexpr tenSquared { 10.0 * 10.0 };
    double constexpr p10_lowLat {1012.0306 - 109.0338 * 10 + 3.6316 * tenSquared};
    const double p72_lowLat {p10_lowLat * exp(-0.147 * (72.0 - 10.0))};

    double constexpr p10_midLat {1012.8186 - 111.5569 * 10.0 + 3.8646 * tenSquared};
    const double p72_midLat {p10_midLat * exp(-0.147 * (72.0 - 10.0))};

    double constexpr p10_highLat {1008.0278 - 113.2494 * 10.0 + 3.9408 * tenSquared};
    const double p72_highLat {p10_lowLat * exp(-0.140 * (72.0 - 10.0))};
}

double GasAttenuationHelpers::convertWaterVaporGM3toHPA(const double &rho_gm3, const double &temp_K) {
    return rho_gm3 * temp_K / 216.7;
}

double GasAttenuationHelpers::calculateImaginaryRefractivity_Oxygen(const double& freq_GHz, 
                const double& dryPressure_hPa, const double& waterVapor_hPa, const double& theta) {
    double lineSum = 0.0;
    double lineStrength, lineShapeFactor;
    for (const auto& coeffRow :
        DataStructures::OXYGEN_COEFFS_TABLE) {
        // Equation #3 ("for oxygen" line) in ITU-R P.676-12
        lineStrength = coeffRow._a1 * 1.0e-7 * dryPressure_hPa * MathHelpers::simpleCube(theta) * std::exp(coeffRow._a2 * (1.0 - theta));
        // Equation #6a ("for oxygen" line) in ITU-R P.676-12
        double lineWidth = coeffRow._a3 * 1.0e-4 * (dryPressure_hPa * std::pow(theta, 0.8 - coeffRow._a4) + 1.1 * waterVapor_hPa * theta);
        double lineWidthSqrd = MathHelpers::simpleSquare(lineWidth);

        // Equation #6b ("for oxygen" line) in ITU-R P.676-12
        lineWidth = std::sqrt(lineWidthSqrd + 2.25e-6);
        lineWidthSqrd = MathHelpers::simpleSquare(lineWidth); // Re-calculate lineWidthSqrd since lineWidth has been updated

        // Equation #7 ("for oxygen" line) in ITU-R P.676-12
        const double correctionFactor = (coeffRow._a5 + coeffRow._a6 * theta) * 1.0e-4 * (dryPressure_hPa + waterVapor_hPa) * std::pow(theta, 0.8);
        const double shapeFactor_part1 = (lineWidth - correctionFactor * (coeffRow._freq_GHz - freq_GHz)) 
                    / (MathHelpers::simpleSquare(coeffRow._freq_GHz - freq_GHz) + lineWidthSqrd);
        const double shapeFactor_part2 = (lineWidth - correctionFactor * (coeffRow._freq_GHz + freq_GHz)) 
                    / (MathHelpers::simpleSquare(coeffRow._freq_GHz + freq_GHz) + lineWidthSqrd);
        lineShapeFactor = (freq_GHz * (shapeFactor_part1 + shapeFactor_part2)) / coeffRow._freq_GHz;
        lineSum += lineStrength * lineShapeFactor;
    }

    // Equation #8 ("for oxygen" line) in ITU-R P.676-12
    const double width_DebyeSpectrum = 5.6e-4 * (dryPressure_hPa + waterVapor_hPa) * std::pow(theta, 0.8);
    const double dryContinuum_part1 = 6.14e-5 / (width_DebyeSpectrum * (1.0 + MathHelpers::simpleSquare(freq_GHz / width_DebyeSpectrum)));
    const double dryContinuum_part2 = 1.4e-12 * dryPressure_hPa * std::pow(theta, 1.5) / (1.0 + 1.9e-5 * std::pow(freq_GHz, 1.5));
    const double dryContinuum = freq_GHz * dryPressure_hPa * MathHelpers::simpleSquare(theta) * (dryContinuum_part1 + dryContinuum_part2);

    double complexRefractivity_oxygen = lineSum + dryContinuum;

    return complexRefractivity_oxygen;
}

double GasAttenuationHelpers::calculateImaginaryRefractivity_Water(const double& freq_GHz, 
                const double& dryPressure_hPa, const double& waterVapor_hPa, const double& theta) {
    double lineSum = 0.0;
    double lineStrength, lineShapeFactor;
    for (const auto& coeffRow :
        DataStructures::WATER_COEFFS_TABLE) {
        // Equation #3 ("for water vapour" line) in ITU-R P.676-12
        lineStrength = coeffRow._b1 * 0.1 * waterVapor_hPa * std::pow(theta, 3.5) * std::exp(coeffRow._b2 * (1.0 - theta));
        // Equation #6a ("for water vapour" line) in ITU-R P.676-12
        double lineWidth = coeffRow._b3 * 1.0e-4 * (dryPressure_hPa * std::pow(theta, coeffRow._b4) + coeffRow._b5 * waterVapor_hPa * std::pow(theta, coeffRow._b6));
        double lineWidthSqrd = MathHelpers::simpleSquare(lineWidth);

        // Equation #6b ("for water vapour" line) in ITU-R P.676-12
        lineWidth = 0.535 * lineWidth + std::sqrt(0.217 * lineWidthSqrd + (2.1316e-12 * coeffRow._freq_GHz) / theta);
        lineWidthSqrd = MathHelpers::simpleSquare(lineWidth);
        
        // Correction factor is set to 0 in Equation #7 ("for water vapor" line) in ITU-R P.676-12, so the following equations have been simplified
        const double shapeFactor_part1 = lineWidth / (MathHelpers::simpleSquare(coeffRow._freq_GHz - freq_GHz) + lineWidthSqrd);
        const double shapeFactor_part2 = lineWidth / (MathHelpers::simpleSquare(coeffRow._freq_GHz + freq_GHz) + lineWidthSqrd);
        lineShapeFactor = (freq_GHz / coeffRow._freq_GHz) * (shapeFactor_part1 + shapeFactor_part2);
        lineSum += lineStrength * lineShapeFactor;
    }

    return lineSum;
}

void GasAttenuationHelpers::setAtmosphericTermsForUsLocation(const double &height_km, 
                double& temp_K, double& totalPressure_hPa, double& waterVapor_hPa, const double &rho0_gm3) {
    if (height_km < 0.0 || height_km > 100.0) {
        std::ostringstream oStrStream; 
        oStrStream << "ERROR: GasAttenuationHelpers::setAtmosphericTermsForUsStandardLocation(): " 
                    << "Cannot determine atmospheric conditions at the given height: " << height_km << " km!";
        throw std::domain_error(oStrStream.str());
    }

    // Convert from height above ground level to geopotential height (pg.1, ITU-R P.835-6)
    const double geopotentialHeight_km = 6356.766 * height_km / (6356.766 + height_km);
    // From page 2, ITU-R P.835-6
    if (height_km <= 86.0) {
        if (geopotentialHeight_km <= 11.0) {
            temp_K = 288.15 - 6.5 * geopotentialHeight_km;
            totalPressure_hPa = 1013.25 * std::pow(288.15 / temp_K, -34.1632 / 6.5);
        }
        else if (geopotentialHeight_km <= 20.0) {
            temp_K = 216.65;
            totalPressure_hPa = 226.3226 * std::exp(-34.1632 * (geopotentialHeight_km - 11.0) / 216.65);
        }
        else if (geopotentialHeight_km <= 32.0) {
            temp_K = 216.65 + geopotentialHeight_km - 20;
            totalPressure_hPa = 54.74980 * std::pow((216.65 / temp_K), 34.1632);
        }
        else if (geopotentialHeight_km <= 47.0) {
            temp_K = 228.65 + 2.8 * (geopotentialHeight_km - 32.0);
            totalPressure_hPa = 8.680422 * std::pow((228.65 / temp_K), (34.1632 / 2.8));
        }
        else if (geopotentialHeight_km <= 51.0) {
            temp_K = 270.65;
            totalPressure_hPa = 1.109106 * std::exp(-34.1632 * (geopotentialHeight_km - 47.0) / 270.65);
        }
        else if (geopotentialHeight_km <= 71.0) {
            temp_K = 270.65 - 2.8 * (geopotentialHeight_km - 51.0);
            totalPressure_hPa = 0.6694167 * std::pow((270.65 / temp_K), (-34.1632 / 2.8));
        }
        else {
            temp_K = 214.65 - 2.0 * (geopotentialHeight_km - 71.0);
            totalPressure_hPa = 0.03956649 * std::pow((214.65 / temp_K), (-34.1632 / 2.0));
        }
    }
    else {
        if (height_km <= 91.0) {
            temp_K = 186.8673;
        }
        else {
            temp_K = 263.1905 - 76.3232 * std::pow((1.0 - MathHelpers::simpleSquare((height_km - 91.0) / 19.9429)), 0.5);
        }
        // Defined in equation #5 (pg.2 of ITU-R P.835-6)
        constexpr double a0 = 95.571899;
        constexpr double a1 = -4.011801;
        constexpr double a2 = 6.424731e-2;
        constexpr double a3 = -4.789660e-4;
        constexpr double a4 = 1.340543e-6;
        totalPressure_hPa = std::exp(a0 + a1 * height_km + a2 * MathHelpers::simpleSquare(height_km) +
                                a3 * MathHelpers::simpleCube(height_km) + a4 * MathHelpers::simpleBiquadrate(height_km));
    }

    constexpr double height0_km = 2.0; // Units of km
    const double rho_gm3 = rho0_gm3 * std::exp(-height_km / height0_km);

    waterVapor_hPa = convertWaterVaporGM3toHPA(rho_gm3, temp_K);
}

void GasAttenuationHelpers::setAtmosphericTermsForUsLowLatitude(const double& height_km, 
                double& temp_K, double& totalPressure_hPa, double& rho_gm3) {
    // Low-latitude annual reference atmosphere
    // Temperature T(k)
    if (height_km >= 0.0 && height_km < 17.0) {
        temp_K = 300.4222 - 6.3533 * height_km + 0.005886 * MathHelpers::simpleSquare(height_km);
    }
    else if (height_km >= 17.0 && height_km < 47.0) {
        temp_K = 194.0 + (height_km - 17.0) * 2.533;
    }
    else if (height_km >= 47.0 && height_km < 52.0) { temp_K = 270.0; }
    else if (height_km >= 52.0 && height_km < 80.0) {
        temp_K = 270.0 - (height_km - 52.0) * 3.0714;
    }
    else if (height_km >= 80.0 && height_km < 100.0) { temp_K = 184.0; }
    // pressure P(hpa)
    if (height_km >= 0.0 && height_km <= 10.0) {
        totalPressure_hPa = 1012.0306 - 109.0338 * height_km + 3.6316 * MathHelpers::simpleSquare(height_km);
    }
    else if (height_km > 10.0 && height_km <= 72) {
        totalPressure_hPa = p10_lowLat * exp(-0.147 * (height_km - 10.0));
    }
    else if (height_km > 72.0 && height_km <= 100.0) {
        totalPressure_hPa = p72_lowLat * exp(-0.165 * (height_km - 72.0));
    }
    // water vapour (g/m^3)
    if (height_km >= 0.0 && height_km <= 15.0) {
        rho_gm3 = 19.6542 * exp((-0.2313 * height_km - 0.1122 * MathHelpers::simpleSquare(height_km) + 0.01351 * pow(height_km, 3.0) - 0.0005923 * pow(height_km, 4.0)));
    }
    else if (height_km > 15.0) { rho_gm3 = 0.0; }
}

void GasAttenuationHelpers::setAtmosphericTermsForUsMidLatitude(const double& height_km, 
                double& temp_K, double& totalPressure_hPa, double& rho_gm3, const Enumerations::Season &season) {
    // Mid-latitude reference atmosphere
    if (season == Enumerations::Season::SummerTime) { // summer mid-latitude
        // Temperature T(k)
        if (height_km < 13.0) {
            temp_K = 294.9838 - 5.2159 * height_km - 0.07109 * pow(height_km, 2);
        }
        else if (height_km < 17.0) { temp_K = 215.15; }
        else if (height_km < 47.0) {
            temp_K = 215.15 * exp((height_km - 17.0) * 0.008128);
        }
        else if (height_km < 53.0) { temp_K = 275.0; }
        else if (height_km < 80.0) {
            temp_K = 275.0 + (1.0 - exp((height_km - 53.0) * 0.06)) * 20.0;
        }
        else { temp_K = 175.0; }
        // pressure P(hPa)
        if (height_km <= 10.0) {
            totalPressure_hPa = 1012.8186 - 111.5569 * height_km + 3.8646 * MathHelpers::simpleSquare(height_km);
        }
        else if (height_km <= 72.0) {
            totalPressure_hPa = p10_midLat * exp(-0.147 * (height_km - 10.0));
        }
        else {
            totalPressure_hPa = p72_midLat * exp(-0.165 * (height_km - 72.0));
        }
        // water vapour (g/m^3)
        if (height_km <= 15.0) {
            rho_gm3 = 14.3542 * exp(-0.4174 * height_km - 0.02290 * MathHelpers::simpleSquare(height_km) + 0.001007 * pow(height_km, 3.0));
        }
        else { rho_gm3 = 0.0; }
    }
    else if (season == Enumerations::Season::WinterTime) { // winter mid-latitude
        const double p10 = 1018.8627 - 124.2954 * 10.0 + 4.8307 * MathHelpers::simpleSquare(10.0);
        const double p72 = p10 * exp(-0.147 * (72.0 - 10.0));
        // Temperature T(k)
        if (height_km < 10.0) {
            temp_K = 272.7241 - 3.6217 * height_km - 0.1759 * pow(height_km, 2);
        }
        else if (height_km < 33.0) { temp_K = 218.0; }
        else if (height_km < 47.0) {
            temp_K = 218.0 + (height_km - 33.0) * 3.3571;
        }
        else if (height_km < 53.0) { temp_K = 265.0; }
        else if (height_km < 80.0) {
            temp_K = 265.0 - (height_km - 53.0) * 2.0370;
        }
        else { temp_K = 210.0; }
        // pressure P(hPa)
        if (height_km <= 10.0) {
            totalPressure_hPa = 1018.8627 - 124.2954 * height_km + 4.8307 * MathHelpers::simpleSquare(height_km);
        }
        else if (height_km <= 72.0) {
            totalPressure_hPa = p10 * exp(-0.147 * (height_km - 10.0));
        }
        else {
            totalPressure_hPa = p72 * exp(-0.155 * (height_km - 72.0));
        }
        // water vapour (g/m^3)
        if (height_km <= 10.0) {
            rho_gm3 = 3.4742 * exp(-0.2697 * height_km - 0.03604 * MathHelpers::simpleSquare(height_km) + 0.0004489 * pow(height_km, 3.0));
        }
        else { rho_gm3 = 0.0; }
    }
}

void GasAttenuationHelpers::setAtmosphericTermsForUsHighLatitude(const double& height_km, 
                double& temp_K, double& totalPressure_hPa, double& rho_gm3, const Enumerations::Season& season) {
    // High latitude reference atmosphere
    if (season == Enumerations::Season::SummerTime) { // summer high-latitude
        // temperature T(k)
        if (height_km < 10.0) {
            temp_K = 286.8374 - 4.7805 * height_km - 0.1402 * MathHelpers::simpleSquare(height_km);
        }
        else if (height_km < 23.0) { temp_K = 225.0; }
        else if (height_km < 48.0) {
            temp_K = 225.0 * exp((height_km - 23.0) * 0.008317);
        }
        else if (height_km < 53.0) { temp_K = 277.0; }
        else if (height_km < 79.0) {
            temp_K = 277.0 - (height_km - 53.0) * 4.0769;
        }
        else { temp_K = 171.0; }
        // pressure (hPa)
        if (height_km <= 10.0) {
            totalPressure_hPa = 1008.0278 - 113.2494 * height_km + 3.9408 * MathHelpers::simpleSquare(height_km);
        }
        else if (height_km <= 72.0) {
            totalPressure_hPa = p10_highLat * exp(-0.140 * (height_km - 10.0));
        }
        else {
            totalPressure_hPa = p72_highLat * exp(-0.165 * (height_km - 72.0));
        }
        // water vapour (g/m^3)
        if (height_km <= 15.0) {
            rho_gm3 = 8.988 * exp(-0.3614 * height_km - 0.005402 * MathHelpers::simpleSquare(height_km) - 0.001955 * pow(height_km, 3.0));
        }
        else { rho_gm3 = 0.0; }
    }
    else if (season == Enumerations::Season::WinterTime) { // winter high-latitude
        // temperature T(k)
        if (height_km < 8.5) {
            temp_K = 257.4345 + 2.3474 * height_km - 1.5479 * MathHelpers::simpleSquare(height_km) + 0.08473 * pow(height_km, 3.0);
        }
        else if (height_km < 30.0) { temp_K = 217.5; }
        else if (height_km < 50.0) { temp_K = 217.5 + (height_km - 30.0) * 2.125; }
        else if (height_km < 54.0) { temp_K = 260.0; }
        else { temp_K = 260.0 - (height_km - 54.0) * 1.667; }
        // pressure (hPa)
        if (height_km <= 10.0) {
            totalPressure_hPa = 1010.8828 - 122.2411 * height_km + 4.554 * MathHelpers::simpleSquare(height_km);
        }
        else if (height_km <= 72.0) {
            const double p10 = 1010.8828 - 122.2411 * 10.0 + 4.554 * MathHelpers::simpleSquare(10.0);
            totalPressure_hPa = p10 * exp(-0.147 * (height_km - 10.0));
        }
        else {
            const double p72 = 1010.8828 - 122.2411 * 72.0 + 4.554 * MathHelpers::simpleSquare(72.0);
            totalPressure_hPa = p72 * exp(-0.150 * (height_km - 72.0));
        }
        //water vapour (g/m^3)
        if (height_km < 10.0) {
            rho_gm3 = 1.2319 * exp(0.07481 * height_km - 0.0981 * MathHelpers::simpleSquare(height_km) + 0.00281 * pow(height_km, 3.0));
        }
        else { rho_gm3 = 0.0; }
    }
}

void GasAttenuationHelpers::setSeasonalAtmosphericTermsForUsLocation(const GeodeticCoord& location, 
                double& temp_K, double& totalPressure_hPa, double& waterVapor_hPa, 
                const Enumerations::Season& season,
                const bool &useAnnualStandardAtmosphere, const double &rho0_gm3) {
    const double height_km = location.height_km;
    const double m_latitude_deg = location.m_latitude_deg;

    if (height_km < 0 || height_km > 100) {
        std::ostringstream oStrStream; 
        oStrStream << "GasAttenuationHelpers::setSeasonalAtmosphericTermsForUsLocation(): The given height is outside of the valid range for heights [0, 100] km: " << height_km << " km!";
        throw std::domain_error(oStrStream.str());
    }
    Enumerations::validateSeason(season);

    double rhoToUpdate_gm3 = 0.0;

    if (useAnnualStandardAtmosphere) {
        // Mean annual global reference atmosphere US Std Atmosphere 1976 
        setAtmosphericTermsForUsLocation(height_km, temp_K, totalPressure_hPa, waterVapor_hPa, rho0_gm3);
        return;
    }
    else if (m_latitude_deg < 22.0) {
        setAtmosphericTermsForUsLowLatitude(height_km, temp_K, totalPressure_hPa, rhoToUpdate_gm3);
        waterVapor_hPa = convertWaterVaporGM3toHPA(rhoToUpdate_gm3, temp_K);
        return;
    }
    else if (m_latitude_deg <= 45.0) {
        setAtmosphericTermsForUsMidLatitude(height_km, temp_K, totalPressure_hPa, rhoToUpdate_gm3, season);
        waterVapor_hPa = convertWaterVaporGM3toHPA(rhoToUpdate_gm3, temp_K);
        return;
    }
    else {
        setAtmosphericTermsForUsHighLatitude(height_km, temp_K, totalPressure_hPa, rhoToUpdate_gm3, season);
        waterVapor_hPa = convertWaterVaporGM3toHPA(rhoToUpdate_gm3, temp_K);
        return;
    }
}

double GasAttenuationHelpers::calculateSpecificWaterAttenuation_dBPerKm(const double& freq_GHz, 
                const double& temp_K, const double& totalPressure_hPa, const double& waterVapor_hPa) {
    // Method in Annex 1 to estimate gaseous attenuation is valid for frequency range 1-1000 GHz (pg.1, ITU-R P.676-12)
    if (freq_GHz < 1.0 || freq_GHz > 1.0e3) {
        std::ostringstream oStrStream;
        oStrStream << "ERROR: GasAttenuation::calculateSpecificOxygenAttenuation_dBPerKm(): The given frequency is outside interval of valid frequencies ([0, 1000] GHz): " << freq_GHz << " GHz!";
        throw std::domain_error(oStrStream.str());
    }

    // Part of Equation #3 from ITU-R P.676-12
    const double theta = 300.0 / temp_K;

    // Equation #3 from ITU-R P.676-12
    const double dryPressure_hPa = totalPressure_hPa - waterVapor_hPa;
    const double complexRefractivity_waterVapor = calculateImaginaryRefractivity_Water(freq_GHz, dryPressure_hPa, waterVapor_hPa, theta);

    const double attenuation_dBPerKm = 0.1820 * freq_GHz * complexRefractivity_waterVapor;

    return attenuation_dBPerKm;
}

double GasAttenuationHelpers::calculateSpecificOxygenAttenuation_dBPerKm(const double& freq_GHz, 
                const double& temp_K, const double& totalPressure_hPa, const double& waterVapor_hPa) {
    // Method in Annex 1 to estimate gaseous attenuation is valid for frequency range 1-1000 GHz (pg.1, ITU-R P.676-12)
    if (freq_GHz < 1.0 || freq_GHz > 1.0e3) {
        std::ostringstream oStrStream;
        oStrStream << "ERROR: GasAttenuation::calculateSpecificWaterAttenuation_dBPerKm(): The given frequency is outside interval of valid frequencies ([0, 1000] GHz): " << freq_GHz << " GHz!";
        throw std::domain_error(oStrStream.str());
    }

    // Part of Equation #3 from ITU-R P.676-12
    const double theta = 300.0 / temp_K;

    // Equation #3 from ITU-R P.676-12
    const double dryPressure_hPa = totalPressure_hPa - waterVapor_hPa;
    const double complexRefractivity_oxygen = calculateImaginaryRefractivity_Oxygen(freq_GHz, dryPressure_hPa, waterVapor_hPa, theta);
    
    const double attenuation_dBPerKm = 0.1820 * freq_GHz * complexRefractivity_oxygen;

    return attenuation_dBPerKm;
}

double GasAttenuationHelpers::calculateSpecificTotalAttenuation_dBPerKm(const double& freq_GHz, 
                const double& temp_K, const double& totalPressure_hPa, const double& waterVapor_hPa) {
    const double oxygenAttenuation_dBPerKm = calculateSpecificOxygenAttenuation_dBPerKm(freq_GHz, temp_K, totalPressure_hPa, waterVapor_hPa);
    const double waterAttenuation_dBPerKm = calculateSpecificWaterAttenuation_dBPerKm(freq_GHz, temp_K, totalPressure_hPa, waterVapor_hPa);

    return oxygenAttenuation_dBPerKm + waterAttenuation_dBPerKm;
}

double GasAttenuationHelpers::calculateRadioRefractiveIndex(const double &totalPressure_hPa, 
                const double &waterVaporPressure_hPa, const double &temp_K) {
    // Equation #6 from ITU-R P.453-14
    const double RADIO_REFRACTIVITY = 77.6 * totalPressure_hPa / temp_K - 5.6 * waterVaporPressure_hPa / temp_K + 3.75e5 * waterVaporPressure_hPa / MathHelpers::simpleSquare(temp_K);
    const double RADIO_REFRACTIVE_IND = 1.0 + RADIO_REFRACTIVITY * 1.0e-6;

    return RADIO_REFRACTIVE_IND;
}


/* ITU-R P.676-12 ANNEX 2 HELPERS START HERE */
double GasAttenuationHelpers::calculateHeightOfWaterComponent_km(const double& freq_GHz, 
                const double& temp_K, const double& rho0_gm3, const double& pressureRatio) {
    // Equation #38
    const double SIGMA_WATER = 1.013 / (1.0 + exp(-8.6 * (pressureRatio - 0.57)));

    double b_factor = 0.0;
    for (uint16_t waterInd = 0; waterInd < DataStructures::WaterFreqList.size(); waterInd++) {
        const double waterA = DataStructures::WaterAList[waterInd];
        const double waterB = DataStructures::WaterBList[waterInd];
        const double waterFreq_GHz = DataStructures::WaterFreqList[waterInd];

        b_factor += (waterA * SIGMA_WATER) / (MathHelpers::simpleSquare(freq_GHz - waterFreq_GHz) + waterB * SIGMA_WATER);
    }

    // Equation #36
    const double A_WATER = 1.9298 - 0.04166 * (temp_K - 273.15) + 0.0517 * rho0_gm3;
    // Equation #37
    const double B_WATER = 1.1674 - 0.00622 * (temp_K - 273.15) + 0.0063 * rho0_gm3;

    // Equation #35(b)
    const double height_waterComp_km = A_WATER + B_WATER * b_factor;

    return height_waterComp_km;
}

double GasAttenuationHelpers::calculateHeightOfOxygenComponent_km(const double& freq_GHz, 
                const double& temp_K, const double& pressureRatio) {
    // Equation #31
    const double T_1_PART1 = 5.1040 / (1.0 + 0.066 * std::pow(pressureRatio, -2.3));
    const double T_1_PART2_SQUARE = (freq_GHz - 59.7) / (2.87 + 12.4 * exp(-7.9 * pressureRatio));
    const double T_1_PART2 = exp(-MathHelpers::simpleSquare(T_1_PART2_SQUARE));
    const double T_1 = T_1_PART1 * T_1_PART2;

    // Equation #32
    double T_2 = 0.0;
    for (uint16_t oxyInd = 0; oxyInd < DataStructures::OxygenConstList.size(); oxyInd++) {
        const double oxyConst = DataStructures::OxygenConstList[oxyInd];
        const double oxyFreq_GHz = DataStructures::OxygenFreqList[oxyInd];
        const double CURRENT_T_2_NUMERATOR = oxyConst * exp(2.12 * pressureRatio);
        const double CURRENT_T_2_DENOMINATOR = MathHelpers::simpleSquare(freq_GHz - oxyFreq_GHz) + 0.025 * exp(2.2 * pressureRatio);
        T_2 += CURRENT_T_2_NUMERATOR / CURRENT_T_2_DENOMINATOR;
    }

    // Equation #33
    const double T_3_PART1 = 0.0114 * freq_GHz / (1.0 + 0.14 * std::pow(pressureRatio, -2.6));
    const double T_3_PART2_NUMERATOR = 15.02 * MathHelpers::simpleSquare(freq_GHz) - 1353.0 * freq_GHz + 5.333e4;
    const double T_3_PART2_DENOMINATOR = MathHelpers::simpleCube(freq_GHz) - 151.3 * MathHelpers::simpleSquare(freq_GHz) + 9629.0 * freq_GHz - 6803.0;
    const double T_3_PART2 = T_3_PART2_NUMERATOR / T_3_PART2_DENOMINATOR;
    const double T_3 = T_3_PART1 * T_3_PART2;

    // Equation #34
    const double A_OXYGEN = 0.7832 + 0.00709 * (temp_K - 273.15);

    // Equation #30
    const double HEIGHT_OXY_PART1 = 6.1 * A_OXYGEN / (1.0 + 0.17 * std::pow(pressureRatio, -1.1));
    double height_oxygenComp_km = HEIGHT_OXY_PART1 * (1.0 + T_1 + T_2 + T_3);

    // Equation #35(a)
    if (freq_GHz < 70.0) {
        height_oxygenComp_km = std::min(height_oxygenComp_km, 10.7 * std::pow(pressureRatio, 0.3));
    }

    return height_oxygenComp_km;
}

double GasAttenuationHelpers::calculateZenithWaterVaporAttenuation_dB(const double &freq_GHz, 
                const double &height_km, const double &integratedWaterVapor_kgm2) {
    if (freq_GHz < 1.0 || freq_GHz > 350.0) {
        std::ostringstream oStrStream;
        oStrStream << "ERROR: GasAttenuation::calculateZenithWaterVaporAttenuation_dB(): The given frequency is outside interval of valid frequencies ([1, 350] GHz): " << freq_GHz << " GHz!";
        throw std::domain_error(oStrStream.str());
    }
    
    // Equation #50
    const double A_PART1 = 0.2048 * exp(-MathHelpers::simpleSquare((freq_GHz - 22.43) / 3.097));
    const double A_PART2 = 0.2326 * exp(-MathHelpers::simpleSquare((freq_GHz - 183.5) / 4.096));
    const double A_PART3 = 0.2073 * exp(-MathHelpers::simpleSquare((freq_GHz - 325.0) / 3.651));
    const double A = A_PART1 + A_PART2 + A_PART3 - 0.1113;

    // Equation #51
    const double B = 8.741e4 * exp(-0.587 * freq_GHz) + 312.2 * std::pow(freq_GHz, -2.38) + 0.723;

    // Equation #52 (restrict height to a value between 0-4 km)
    double calcHeight_km = 0.0;
    if (height_km >= 0.0 && height_km <= 4.0) {
        calcHeight_km = height_km;
    }
    else {
        calcHeight_km = 4.0;
    }

    // Equation #53-54
    const double RHO_VREF_GM3 = integratedWaterVapor_kgm2 / 2.38;
    const double TEMP_REF_C = 14.0 * std::log(0.22 * integratedWaterVapor_kgm2 / 2.38) + 3.0;
    const double TEMP_REF_K = TEMP_REF_C + 273.15;
    const double WATER_PRESSURE_REF_HPA = convertWaterVaporGM3toHPA(RHO_VREF_GM3, TEMP_REF_K);
    const double FREQ_REF_GHZ = 20.6;
    const double DRY_PRESSURE_REF_HPA = 845.0;
    const double TOTAL_PRESSURE_REF_HPA = DRY_PRESSURE_REF_HPA + WATER_PRESSURE_REF_HPA;

    // Equation #49
    const double GAMMA_WATER = calculateSpecificWaterAttenuation_dBPerKm(freq_GHz, TEMP_REF_K, TOTAL_PRESSURE_REF_HPA, WATER_PRESSURE_REF_HPA);
    const double GAMMA_REF_WATER = calculateSpecificWaterAttenuation_dBPerKm(FREQ_REF_GHZ, TEMP_REF_K, TOTAL_PRESSURE_REF_HPA, WATER_PRESSURE_REF_HPA);

    const double WATER_VAPOR_ATTEN_DB = 0.0176 * integratedWaterVapor_kgm2 * GAMMA_WATER / GAMMA_REF_WATER;

    if (freq_GHz < 20.0) {
        return WATER_VAPOR_ATTEN_DB;
    }
    else {
        const double FACTOR = A * std::pow(calcHeight_km, B) + 1.0;
        return WATER_VAPOR_ATTEN_DB * FACTOR;
    }
}