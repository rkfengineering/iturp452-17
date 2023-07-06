#include "MathHelpers.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

double MathHelpers::rootRaisedCosine(const double& t, const double& alpha)
{
	double z;

	z = alpha * (std::cos(M_PI * (t + 0.25)) * sinc(alpha * t + 0.25)
		+ std::cos(M_PI * (t - 0.25)) * sinc(alpha * t - 0.25))
		+ (1.0 - alpha) * sinc(t * (1.0 - alpha));

	return(z);
};

double MathHelpers::interpolate1D(const double& startValue, const double& endValue, const double& endValueWeight) {
	return startValue + endValueWeight * (endValue - startValue);
}

/// <summary>
/// Calculate the K-weight described in Section 2 "Bi-cubic interpolation" from ITU-R P.1144-6.
/// This is the first step to calculate the actual bi-cubic interpolated value.
/// </summary>
/// <param name="delta">Input to the K-weight function</param>
/// <returns>K-weight for use in calculating bi-cubic interpolation</returns>
double MathHelpers::calculateBicubicInterpolationWeight(const double& delta) {
	const double A = -0.5;
	const double ABS_DELTA = std::abs(delta);

	if (ABS_DELTA <= 1.0) {
		return (A + 2.0) * simpleCube(ABS_DELTA) - (A + 3.0) * simpleSquare(ABS_DELTA) + 1.0;
	}
	else if (ABS_DELTA <= 2.0) {
		return A * simpleCube(ABS_DELTA) - 5.0 * A * simpleSquare(ABS_DELTA) + 8.0 * A * ABS_DELTA - 4.0 * A;
	}
	else {
		return 0.0;
	}
}

double MathHelpers::calculateBicubicInterpolation(const std::vector<std::vector<double>>& gridValueMatrix, const double& rowWeight, const double& columnWeight) {
	// Given: Values on a 4x4 grid
	const uint16_t GRID_AXIS_LENGTH = 4;

	if (gridValueMatrix.size() != GRID_AXIS_LENGTH || gridValueMatrix[0].size() != GRID_AXIS_LENGTH) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: MathHelpers::calculateBicubicInterpolation(): " 
					<< "Provided grid does not have the correct dimensions of 4x4: " 
					<< gridValueMatrix.size() << "x" << gridValueMatrix[0].size() << "!";
		throw std::domain_error(oStrStream.str());
	}

	// Step 1: For each row compute the interpolated value
	std::vector<double> colScaleFactorList;
	colScaleFactorList.reserve(GRID_AXIS_LENGTH);
	double cubicInterpWeight;
	for (uint16_t colInd = 0; colInd < GRID_AXIS_LENGTH; ++colInd) {
		cubicInterpWeight = MathHelpers::calculateBicubicInterpolationWeight(colInd - (columnWeight + 1.0));
		colScaleFactorList.push_back(cubicInterpWeight);
	}
	std::vector<double> interpRowList;
	interpRowList.reserve(GRID_AXIS_LENGTH);
	double interpRowValue;
	for (uint16_t rowInd = 0; rowInd < GRID_AXIS_LENGTH; rowInd++) {
		interpRowValue = 0.0; // Initialize to zero
		for (uint16_t colInd = 0; colInd < GRID_AXIS_LENGTH; colInd++) {
			interpRowValue += gridValueMatrix[rowInd][colInd] * colScaleFactorList[colInd];
		}
		interpRowList.push_back(interpRowValue);
	}

	// Step 2: Interpolate the one-dimensional interpolations in the same manner as the row interpolations
	double interpResult = 0.0;
	for (uint16_t rowInd = 0; rowInd < GRID_AXIS_LENGTH; rowInd++) {
		cubicInterpWeight = MathHelpers::calculateBicubicInterpolationWeight(rowInd - (rowWeight + 1.0));
		interpResult += interpRowList[rowInd] * cubicInterpWeight;
	}

	return interpResult;
}

/// <summary>
/// Assume that xList and yList are associated element-by-element. By using the placement of the inputX amongst xList, determine the associated yValue by interpolating from the yList.
/// </summary>
/// <param name="inputX">x-value of the desired y-value</param>
/// <param name="xList">List of x-values to extract weight from</param>
/// <param name="yList">List of y-values to interpolate within</param>
/// <returns>Interpolated y-value associated with the given x-value</returns>
double MathHelpers::interpolateAlongXForY(double inputX, std::vector<double> xList, std::vector<double> yList) {
	auto xIter = std::adjacent_find(xList.begin(), xList.end(), std::greater_equal<double>());
	if (xIter != xList.end()) {
		throw std::logic_error("ERROR: MathHelpers::interpolate2D(): " 
					"Input list of x-values are not strictly sorted in ascending order, "
					"therefore this interpolation method will not work!");
	}

	// Don't perform extrapolation if the input value is greater than the largest value in the vector
	if (inputX >= xList.back()) {
		return yList.back();
	}
	// Same logic again if the input value is less than the smallest value in the vector
	if (inputX <= xList.front()) {
		return yList.front();
	}

	const auto GT_OR_EQUAL_ITER = std::lower_bound(xList.begin(), xList.end(), inputX); // First element in the vector where: value >= inputValue
	const uint64_t UPPER_VALUE_POSITION = (GT_OR_EQUAL_ITER - xList.begin());

	const double UPPER_X_VALUE = xList[UPPER_VALUE_POSITION];
	const double LOWER_X_VALUE = xList[UPPER_VALUE_POSITION - 1];
	const double UPPER_Y_VALUE = yList[UPPER_VALUE_POSITION];
	const double LOWER_Y_VALUE = yList[UPPER_VALUE_POSITION - 1];

	const double UPPER_VALUE_WEIGHT = (inputX - LOWER_X_VALUE) / (UPPER_X_VALUE - LOWER_X_VALUE);
	const double RESULT = interpolate1D(LOWER_Y_VALUE, UPPER_Y_VALUE, UPPER_VALUE_WEIGHT);

	return RESULT;
}

double MathHelpers::unwrapValueAroundAxis(const double& value, const double& minValue, const double& maxValue) {
	// For example, if the axis bounds are 10 --> 100, then a value of 5 will result in floor(-5/90) = -1 wraps
	// -1 wraps will change that value from 5 to 95 (5 - (90 * -1) = 95)
	const double AXIS_SIZE = maxValue - minValue;
	// This is the number of (positive or negative) wraps occurring
	const double NUM_WRAPS = std::floor((value - minValue) / AXIS_SIZE);
	// Remove the number of wraps
	const double UNWRAPPED_VALUE = value - AXIS_SIZE * NUM_WRAPS;
	return UNWRAPPED_VALUE;
}

double MathHelpers::clampValueWithinAxis(const double& value, const double &maxValue, const double &minValue) {
	const double CLAMP_TO_AXIS_MAX = std::min(maxValue, value);
	const double CLAMP_TO_AXIS_MIN = std::max(minValue, CLAMP_TO_AXIS_MAX);
	return CLAMP_TO_AXIS_MIN;
}

double MathHelpers::mirrorValueAroundAxis(const double &value, const double& maxValue, const double& minValue) {
	// This is the number of (positive or negative) wraps occurring
	const double NUM_WRAPS_FLOAT = std::floor(value / maxValue);
	const int32_t NUM_WRAPS = static_cast<int32_t>(NUM_WRAPS_FLOAT);
	double unwrappedValue = 0.0;
	// If it wraps around the axis an even number of times, there's no need to mirror it
	if (NUM_WRAPS % 2 == 0) {
		// Remove the number of wraps
		unwrappedValue = value - maxValue * NUM_WRAPS_FLOAT;
	}
	// If it wraps around the axis an odd number of times, it must be inverted within the axis in order to mirror it
	else {
		unwrappedValue = maxValue * NUM_WRAPS_FLOAT - value;
	}
	return unwrappedValue;
}

double MathHelpers::getAverage(const std::vector<double> &valueList)
{
	if (valueList.empty()) {
		return 0.0;
	}

	const double AVERAGE = std::accumulate(valueList.begin(), valueList.end(), 0.0) / valueList.size();
	return AVERAGE;
}

double MathHelpers::getProductAverage(const std::vector<double> &firstList, const std::vector<double> &secondList)
{
	if (firstList.size() != secondList.size()) {
		std::ostringstream oStrStream;
		oStrStream << "ERROR: MathHelpers::getProductAverage(): Provided lists are not of equal length: " 
					<< firstList.size() << " != " << secondList.size() << ". Cannot perform elemnt-wise multiplication!";
		throw std::domain_error(oStrStream.str());
	}

	std::vector<double> productVector;
	productVector.reserve(firstList.size());
	// Element-wise multiple of first and second lists
	std::transform(firstList.begin(), firstList.end(),
		secondList.begin(), std::back_inserter(productVector),
		std::multiplies<double>());

	// Calculate average of the product vector
	return getAverage(productVector);
}

double MathHelpers::getLeastSquaresBestFitSlope(const std::vector<double>& xAxis, const std::vector<double>& yAxis)
{
	const double NUMERATOR = MathHelpers::getProductAverage(xAxis, yAxis) - MathHelpers::getAverage(xAxis) * MathHelpers::getAverage(yAxis);
	const double DENOMINATOR = MathHelpers::getProductAverage(xAxis, xAxis) - MathHelpers::getAverage(xAxis) * MathHelpers::getAverage(xAxis);
	const double MIN_DENOMINATOR = 1.0e-20;
	const double RESULT = NUMERATOR / std::max(MIN_DENOMINATOR, DENOMINATOR);

	return RESULT;
}

double MathHelpers::getYInterceptOfLeastSquaresBestFit(const std::vector<double> &xAxis, const std::vector<double> &yAxis)
{
	const double X_AXIS_DOUBLE_PRODUCT_AVERAGE = MathHelpers::getProductAverage(xAxis, xAxis);
	const double X_AXIS_AVERAGE = MathHelpers::getAverage(xAxis);

	const double NUMERATOR = X_AXIS_DOUBLE_PRODUCT_AVERAGE * MathHelpers::getAverage(yAxis) - MathHelpers::getProductAverage(xAxis, yAxis) * X_AXIS_AVERAGE;
	const double DENOMINATOR = X_AXIS_DOUBLE_PRODUCT_AVERAGE - X_AXIS_AVERAGE * X_AXIS_AVERAGE;
	const double MIN_DENOMINATOR = 1.0e-20;

	const double RESULT = NUMERATOR / std::max(MIN_DENOMINATOR, DENOMINATOR);
	return RESULT;
}

int MathHelpers::getGreatestCommonDivisor(const int &firstValue, const int &secondValue) {
	// EXAMPLE: firstValue = -25, secondValue = 5
	// STEP 0 --> firstValue = 5, secondValue = -25 % 5 = 0
	// STEP 1 --> secondValue == 0, returns firstValue = 5
	return secondValue == 0 ? firstValue : getGreatestCommonDivisor(secondValue, firstValue % secondValue);
}


template<typename T>
bool MathHelpers::areSameOrBothNaN(const T& firstValue, const T& secondValue) {
	// If both are defined and equal, consider them equal
	if (firstValue == secondValue) {
		return true;
	}
	// If both are NaN, consider them equal
	return (std::isnan(firstValue) && std::isnan(secondValue));
}

double MathHelpers::Qnorm(const double &x)
{
	double z = 0.50 * x * std::sqrt(2.0);
	z = std::erf(z);
	z = 0.50 - 0.50 * z;
	return z;
}

double MathHelpers::invQnorm(const double& prob)
{
	const double CLAMPED_PROB = prob;
	// We need to clamp probability between 0 and 1.
	// So any value greater than 1 will become 1, and any value lower than 0 will become 0.
	// Now, if the probability is exactly 0 or 1, we already know the answer
	if (CLAMPED_PROB >= 1.0)
		return 0.0;
	else if (CLAMPED_PROB <= 0.0)
		return 1.0;

	double actualProb = 1.0 - 2.0 * prob;
	double w = 1.0 - actualProb * actualProb;
	w = -1.0 * std::log(w);
	double x = 0.0;
	if (w < 5.0)
	{
		w -= 2.5;
		x = 2.81022636e-08;
		x = 3.43273939e-07 + x * w;
		x = -3.5233877e-06 + x * w;
		x = -4.39150654e-06 + x * w;
		x = 0.00021858087 + x * w;
		x = -0.00125372503 + x * w;
		x = -0.00417768164 + x * w;
		x = 0.246640727 + x * w;
		x = 1.50140941 + x * w;
	}
	else
	{
		w = std::sqrt(w) - 3.0;
		x = -0.000200214257;
		x = 0.000100950558 + x * w;
		x = 0.00134934322 + x * w;
		x = -0.00367342844 + x * w;
		x = 0.00573950773 + x * w;
		x = -0.0076224613 + x * w;
		x = 0.00943887047 + x * w;
		x = 1.00167406 + x * w;
		x = 2.83297682 + x * w;
	}

	double z = actualProb * x;
	z = std::copysign(z, actualProb);
	z = std::sqrt(2.0) * z;

	return z;
}

double MathHelpers::invBivarQnorm(const double &prob, const double &minXValue, const double &correlation)
{
	double clampedCorrelation = correlation;
	if (correlation > 1.0) clampedCorrelation = 1.0;
	if (correlation < -1.0) clampedCorrelation = -1.0;
	double clampedProb = prob;
	if (prob > 1.0) clampedProb = 1.0;
	if (prob < 0.0) clampedProb = 0.0;
	
	const uint32_t MAX_NUM_RUNS = 100;
	uint32_t runInd = 0;
	double invResult = 0.0;
	double H = 20.0;
	double L = -20.0;
	double bivarQnormResult = bivarQnorm(minXValue, invResult, clampedCorrelation);
	double diffFromDesiredProb = std::abs(bivarQnormResult - clampedProb);

	if (clampedProb > Qnorm(minXValue)) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	// Recursively calculate value until the desired probability is reached, or until we've gotten as close as possible
	while (diffFromDesiredProb > 1.0e-12 && runInd < MAX_NUM_RUNS)
	{
		if (bivarQnormResult > clampedProb) {
			L = invResult;
			invResult = (invResult + H) / 2.0;
		}
		else {
			H = invResult;
			invResult = (invResult + L) / 2.0;
		}
		bivarQnormResult = bivarQnorm(minXValue, invResult, clampedCorrelation);
		diffFromDesiredProb = std::abs(bivarQnormResult - clampedProb);
		runInd++;
	}

	return invResult;
}

double MathHelpers::bivarQnorm(const double& a, const double& b, const double& r) {
	if ((a > 10.0) || (b > 10.0)) {
		return 0.0;
	}

	if ((a < -10.0)) {
		return Qnorm(b);
	}

	if ((b < -10.0)) {
		return Qnorm(a);
	}

	if (r >= 1.0) {
		return Qnorm(std::max(a, b));
	}

	if (r <= -1.0) {
		if (a + b >= 0) {
			return 0.0;
		}
		else {
			return 1.0 - (Qnorm(-a) + Qnorm(-b));
		}
	}

	double h1 = -1.0 * a;
	double h2 = -1.0 * b;
	double h12 = (h1 * h1 + h2 * h2) / 2.0;

	double h3, h5, h6, h7, h8;
	double A;
	double ab;
	double LH;
	double r1, r2, r3;
	double rr;
	double x[] = {0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992};
	double w[] = {0.018854042, 0.038088059, 0.0452707394, 0.038088059,
				  0.018854042};
	double rAbs = std::abs(r);

	double B = 0.0; // This initial value will get reset when necessary
	if (rAbs >= 0.7) {
		r2 = 1.0 - r * r;
		r3 = std::sqrt(r2);
		if (r < 0.0) {
			h2 = -h2;
		}
		h3 = h1 * h2;
		h7 = std::exp(-0.50 * h3);
		
		// By this point, we know that rAbs < 1
		h6 = std::abs(h1 - h2);
		h5 = 0.5 * h6 * h6;
		h6 = h6 / r3;
		A = 0.5 - h3 / 8.0;
		ab = 3.0 - 2.0 * A * h5;
		LH = 0.13298076 * h6 * ab * (Qnorm(h6)) -
				std::exp(-1.0 * h5 / r2) * (ab + A * r2) * 0.053051647;
		for (uint16_t j = 0; j < 5; j++)
		{
			r1 = r3 * x[j];
			rr = r1 * r1;
			r2 = std::sqrt(1 - rr);
			if (h7 == 0)
				h8 = 0;
			else
				h8 = std::exp(-1.0 * h3 / (1.0 + r2)) / (r2 * h7);
			LH = LH - w[j] * std::exp(-1.0 * h5 / rr) * (h8 - 1.0 - A * rr);
		}
		B = h1;
		if (h2 < h1) {
			B = h2;
		}
		B = LH * r3 * h7 + 1.0 - Qnorm(B);
		if (r < 0.0) {
			B = 1.0 - Qnorm(h1) - B;
		}

		return B;
	}
	
	// If rAbs < 0.7, then follow this logic
	h3 = h1 * h2;
	LH = 0;
	if (r != 0)
	{
		for (uint16_t j = 0; j < 5; j++)
		{
			r1 = r * x[j];
			r2 = 1.0 - r1 * r1;
			LH = LH + w[j] * std::exp((r1 * h3 - h12) / r2) / std::sqrt(r2);
		}
	}

	B = (1.0 - Qnorm(h1)) * (1.0 - Qnorm(h2)) + r * LH;

	return B;
}