#ifndef MATH_HELPERS_H
#define MATH_HELPERS_H

#define _USE_MATH_DEFINES // to allow M_PI definition

#include <math.h>
#include <vector>

namespace MathHelpers {
	/** Round a float to nearest integer value.
	 *
	 * @param x The value to round.
	 * @return The rounded value.
	 */
	template<typename T>
	T round(T x) {
		return std::floor(x + 0.5);
	}

	/** Convert an angle from degrees to radians.
	 * @param deg The angle in degrees
	 * @return The angle in radians.
	 */
	template<typename T>
	T deg2rad(T deg) {
		return M_PI / 180.0 * deg;
	}
	/** Convert an angle from radians to degrees.
	 * @param rad The angle in radians
	 * @return The angle in degrees.
	 */
	template<typename T>
	T rad2deg(T rad) {
		return 180.0 / M_PI * rad;
	}

	/** Shortcut for computing squares.
	 * @param val The value to square.
	 * @return The square of val.
	 */
	template<typename T>
	T simpleSquare(T val) {
		return val * val;
	}

	/** Shortcut for computing cubes.
	 * @param val The value to cube.
	 * @return The cube of val.
	 */
	template<typename T>
	T simpleCube(T val) {
		return val * val * val;
	}

	/** Shortcut for computing x^4.
	 * @param val The value to x^4.
	 * @return The cube of val.
	 */
	template<typename T>
	T simpleBiquadrate(T val) {
		return val * val * val * val;
	}

	/** Shortcut for computing sinc.
	 * @param x The value to compute sinc(x) = sin(pi * x) / (pi * x).
	 * @return The sinc of @c x.
	 */
	template <typename T>
	T sinc(T x) {
		static const double eps = 1e-6;
		if (x < eps && x > -eps) {
			return 1.0 - MathHelpers::simpleSquare(M_PI * x) / 6.0;
		}
		else {
			return std::sin(M_PI * x) / (M_PI * x);
		}
	}

	/** Root-raised cosine pulse.
	 * @param t time of root-raised cosine sample.
	 * @param alpha spectral rolloff.
	 * @return The sample of root-raised cosine pulse.
	 */
	double rootRaisedCosine(const double& t, const double& alpha);

	/// <summary>
	/// Interpolate between two points, given a weight
	/// </summary>
	/// <param name="startValue">Start point for interpolation</param>
	/// <param name="endValue">End point for interpolation</param>
	/// <param name="endValueWeight">Weight applied to the end value (weight = 1 --> result == endValue, weight = 0 --> result == startValue)</param>
	/// <returns>Interpolated value between start and end values</returns>
	double interpolate1D(const double& startValue, const double& endValue, const double& endValueWeight);

	/// <summary>
	/// Calculate the K-weight described in Section 2 "Bi-cubic interpolation" from ITU-R P.1144-6.
	/// This is the first step to calculate the actual bi-cubic interpolated value.
	/// </summary>
	/// <param name="delta">Input to the K-weight function</param>
	/// <returns>K-weight for use in calculating bi-cubic interpolation</returns>
	double calculateBicubicInterpolationWeight(const double& delta);

	double calculateBicubicInterpolation(const std::vector<std::vector<double>>& gridValueMatrix, const double& rowWeight, const double& columnWeight);

	/// <summary>
	/// Assume that xList and yList are associated element-by-element. By using the placement of the inputX amongst xList, determine the associated yValue by interpolating from the yList.
	/// </summary>
	/// <param name="inputX">x-value of the desired y-value</param>
	/// <param name="xList">List of x-values to extract weight from</param>
	/// <param name="yList">List of y-values to interpolate within</param>
	/// <returns>Interpolated y-value associated with the given x-value</returns>
	double interpolateAlongXForY(double inputX, std::vector<double> xList, std::vector<double> yList);

	double unwrapValueAroundAxis(const double& value, const double& minValue, const double& maxValue);

	double clampValueWithinAxis(const double& value, const double &maxValue, const double &minValue = 0.0);

	/** Wrap a value to a particular size by mirroring the object space onto
	 * the image space.
	 * @param size The exclusive maximum limit.
	 * @param value The value to be limited.
	 * @return The value limited to the range [0, @a size) by mirroring.
	 */
	double mirrorValueAroundAxis(const double &value, const double& maxValue, const double& minValue = 0.0);

	/** Returns the average of a a vector
	 *
	 * @param arr[] 1-d vector of values to be averaged
	 * @return The average value of the values in arr
	 */
	double getAverage(const std::vector<double> &valueList);
	
	/** Returns the average of the element-by-element product of two equally-sized 1-d vectors
	*
	* @param arr1[] First 1-d vector of values
	* @param arr2[] Second 1-d vector of values
	* @return The product-average value of the values in the two vectors
	*/
	double getProductAverage(const std::vector<double> &firstList, const std::vector<double> &secondList);

	/** Returns the linear-best-fit least-squares slope for the data (xs,ys)
	*
	* @param xs a 1-d vector of x-values
	* @param ys a 1-d vector of y-values
	* @return The linear-best-fit slope
	*/
	double getLeastSquaresBestFitSlope(const std::vector<double>& xAxis, const std::vector<double>& yAxis);

	/** Returns the linear-best-fit least-squares y-intercept for the data (xs,ys)
	*
	* @param xs a 1-d vector of x-values
	* @param ys a 1-d vector of y-values
	* @return The linear-best-fit y-intercept
	*/
	double getYInterceptOfLeastSquaresBestFit(const std::vector<double> &xAxis, const std::vector<double> &yAxis);

	struct NeighborIntegerPair {
		/// <summary>
		/// Store a double value in between two bounding integer-value doubles, calculating a weight factor for each bounding double
		/// </summary>
		/// <param name="value">Value to be stored</param>
		NeighborIntegerPair(const double &value) {
			lowPoint = std::floor(value);
			highPoint = std::ceil(value);
			weightFactor = value - lowPoint;
		}
		// First integer number below the given value
		double lowPoint;
		// First integer number above the given value
		double highPoint;
		// Weight of highPoint (weightFactor = 1 --> value == highPoint, weightFactor = 0, value == lowPoint)
		double weightFactor;
	};

	/**
	 * Compute GCD of two integers
	 */
	int getGreatestCommonDivisor(const int &firstValue, const int &secondValue);

	/// <summary>
	/// Determine if two numbers are equal, additionally accounting for the possibility that both values are NaN (in which case they are equal too)
	/// </summary>
	/// <typeparam name="T">Data type of the numbers being compared</typeparam>
	/// <param name="firstValue">First value to compare</param>
	/// <param name="secondValue">Second value to compare</param>
	/// <returns>Boolean indicating if the two numbers are equal, or both NaN</returns>
	template<typename T>
	bool areSameOrBothNaN(const T& firstValue, const T& secondValue);

	/** Q function (probability a standard normal variable exceeds a value)
	* Q(x)=0.5-0.5*erf(x/sqrt(2))
	* @param x the value to be exceeded
	* @return The probability a standard normal variable exceeds that value
	*/
	double Qnorm(const double &x);
	
	/** Inverse Q function (the value required so that the probability a
	 * standard normal variable exceeds it is a given value)
	 * Qinverse(p)=sqrt(2)*erfinverse(2*p-1)
	 * From Mike Giles' "Approximating the erfinv function"
	 *
	 * @param p desired probability
	 * @return The required value to achieve the probability
	 */
	double invQnorm(const double& prob);

	/** Bivariate Q function (returns the probability that a bivariate-standard-normal random variable with
	 * correlation r will lie in the region [a,infinity]x[b,infinity]
	 * From Graeme West's "Better Approximations to Cumulative Normal Functions"
	 *
	 * @param a x-value lower bound
	 * @param b y-value lower bound
	 * @param r correlation (abs(r)<=1)
	 * @return The probability that the random variable will be in [a,infinity]x[b,infinity]
	 */
	double bivarQnorm(const double& a, const double& b, const double& r);

	/** Inverse bivariate Q function (returns the value b such that a bivariate-standard-normal random variable with
	 * correlation r will lie in the region [a,infinity]x[b,infinity] with probability p.
	 * Based on a binary search algorithm. Returns NaN if the input probability is unobtianable
	 * for any value of b given the value of a (i.e. if p>Qnorm(a)).
	 *
	 * @param p Desired probability.
	 * @param a x-value lower bound
	 * @param r correlation (abs(r)<=1)
	 * @return The requisite x-value lower bound.
	 */
	double invBivarQnorm(const double &prob, const double &minXValue, const double &correlation);
} // end namespace MathHelpers
#endif /* MATH_HELPERS_H */