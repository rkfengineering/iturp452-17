#ifndef INV_CUM_NORM_H
#define INV_CUM_NORM_H

//TODO figure out where to put this. maybe add to mathhelpers

//From ITU-R P.452-17 Annex 1 Attachment 3

/// @brief Approximation to inverse cumulative normal distribution function for x<0.5
/// @param x threshold probability (fraction)
/// @return Value from normal distribution (mu = 0, sigma = 1)
double inv_cum_norm(double x);

#endif /* INV_CUM_NORM_H */