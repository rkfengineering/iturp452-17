
#ifndef ENUMS_H
#define ENUMS_H

#include <sstream>
#include <iostream>

namespace Enumerations {

	/// <summary>
	/// Specifies signal polarization type.
	/// First value is set to 1 to require that this value is set intentionally (default value 0 has no meaning).
	/// </summary>
	enum PolarizationType {
		HorizontalPolarized = 1,
		VerticalPolarized = 2,
		CircularPolarized = 3,
	};

	/// <summary>
	/// Specifies season for atmospheric modelling.
	/// First value is set to 1 to require that this value is set intentionally (default value 0 has no meaning).
	/// </summary>
	enum Season {
		SummerTime = 1,
		WinterTime = 2,
	};

	/// <summary>
	/// Determine the signal's angle given the polarization of an antenna.
	/// </summary>
	/// <param name="polarType">Polarization type of the antenna</param>
	/// <returns>Signal's angle (deg)</returns>
	const double determineAngleOfPolarization_deg(const PolarizationType& polarType);

	/// <summary>
	/// Throw an error if a season is not valid.
	/// </summary>
	/// <param name="season">Season type to validate</param>
	/// <returns>Nothing. Throws if the season type is not valid</returns>
	void validateSeason(const Season& season);
}

#endif /* ENUMS_H */
