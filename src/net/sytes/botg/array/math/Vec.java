package net.sytes.botg.array.math;

import java.util.Arrays;

import net.sytes.botg.array.ArUtils;
import net.sytes.botg.array.ConvertArray;
import net.sytes.botg.array.SearchArray;

public class Vec {

	public enum DownsamplingAlgorithm {
		BRUTE_FORCE, MAX, MEAN, LARGEST_TRIANGLE_THREE_BUCKETS, LARGEST_TRIANGLE_ONE_BUCKET, LONGEST_LINE_BUCKET
	}
	
	// Suppress default constructor for noninstantiability
	private Vec() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * apply a highpass on {@code x} defined by {@code dt} and {@code fc}
	 * @param x
	 * @param dt equidistant time period between samples
	 * @param f_c cutoff frequency
	 * @return
	 */
	public static double[] highpass(double[] x, double dt, double f_c) {
		int L = x.length;
		double[] x_f = new double[L];
		double alpha = 1 / (1 + 2 * Math.PI * f_c * dt);
		x_f[0] = x[0];
		for (int i = 1; i < L; i++) {
			x_f[i] = alpha * x_f[i - 1] + alpha * (x[i] - x[i - 1]);
		}
		return x_f;
	}
	
	/**
	 * apply a lowpass on {@code x} defined by {@code dt} and {@code fc}
	 * @param x
	 * @param dt equidistant time period between samples
	 * @param f_c cutoff frequency
	 * @return
	 */
	public static double[] lowpass(double[] x, double dt, double f_c) {
		int L = x.length;
		double[] x_f = new double[L];
		double alpha = dt / (dt + 1 / (2 * Math.PI * f_c));
		x_f[0] = alpha * x[0];
		for (int i = 1; i < L; i++) {
			x_f[i] = alpha * x[i] + (1 - alpha) * x_f[i - 1];
		}
		return x_f;
	}
	
	/**
	 * TODO not implemented
	 * 1D linear interpolation
	 * @param x original x-value vector
	 * @param y original y-value vector
	 * @param xi new x-value vector to interpolate for
	 * @return yi interpolated y-values as double[] 
	 */
	public static double[] interp1Lin(double[] x, double[] y, double[] xi) {
		// TODO finish impl
		ArUtils.checkForNull(x);
		ArUtils.checkForEmpty(y);
		ArUtils.checkForNull(xi);
		ArUtils.checkForEqualDimensions(x, y);
		int n = x.length;
		int ni = xi.length;
		double[] yi = new double[ni];
		
		int i = 0;
		int j = 0;
		
		double dx = 0.0;
		double dy = 0.0;
		double dxi = 0.0;
		
		// check if xi starts lower than x, then interpolate accordingly
		while(xi[i] <= x[j]) {
			dx = x[j + 1] - x[j];
			dy = y[j + 1] - y[j];
			dxi = x[j] - xi[i];
			yi[i] = y[i] - dy * dxi / dx;
			++i;
		}
		
		++j;
		
		// go through rest
		while (i < n) {
			
			// find closest lower x
			
			// find closest upper x
			
			++i;
		}
		
		return yi;
	}
	
	/**
	 * removes {@code start} samples from the beginning of the {@code ar} and {@code end} samples from the end of {@code ar} 
	 * <br>and returns the result as new double[] array 
	 * @param start
	 * @param end
	 * @return
	 */
	public static double[] trim(double[] ar, int start, int end) {
		ArUtils.checkForGreaterEqualZero(new int[]{start, end});
		int n = ar.length;
		if (start + end > n) {
			throw new IllegalArgumentException("sum of start and end must be smaller than length of ar");
		}
		return ArUtils.subArray(ar, start, n - end); 
	}
	
	/**
	 * normalization of vector {@code x} into range [0, 1]
	 * @param x
	 * @return
	 */
	public static double[] normalize(double[] x) {
		double xMax = Vec2Scalar.max(x);
		double xMin = Vec2Scalar.min(x);
		int n = x.length;
		double[] xn = new double[n];
		for (int i = 0; i < n; i++) {
			xn[i] = (x[i] - xMin) / (xMax - xMin);
		}
		return xn;
	}
	
	/**
	 * normalization of vector {@code x} into range [-1, 1]
	 * @param x
	 * @return
	 */
	public static double[] normalize2(double[] x) {
		double xMax = Vec2Scalar.max(x);
		double xMin = Vec2Scalar.min(x);
		int n = x.length;
		double[] xn = new double[n];
		for (int i = 0; i < n; i++) {
			xn[i] = 2 * ((x[i] - xMin) / (xMax - xMin) - 0.5);
		}
		return xn;
	}
	
	/**
	 * normalization of vector {@code x} with zscore
	 * @param x
	 * @return
	 */
	public static double[] zscore(double[] x) {
		double mean = Vec2Scalar.mean(x);
		double sigma = Vec2Scalar.variance(x);
		int n = x.length;
		double[] z = new double[n];
		for (int i = 0; i < n; i++) {
			z[i] = (x[i] - mean) / sigma;
		}
		return z;
	}
	
	/**
	 * downsamples an array using the specified {@code algorithm}
	 * <br><a href="https://skemman.is/bitstream/1946/15343/3/SS_MSthesis.pdf">Link</a>
	 * @param ar
	 * @param threshold
	 * @param algorithm Options -> BRUTE_FORCE | MAX | MEAN | LARGEST_TRIANGLE_THREE_BUCKETS | LARGEST_TRIANGLE_ONE_BUCKET | LONGEST_LINE_BUCKET
	 * @return double[]
	 */
	public static double[] downsample(double[] ar, int m, DownsamplingAlgorithm algorithm) {		
		switch(algorithm) {
			case BRUTE_FORCE:
				return downsampleBruteForce(ar, ar.length / m);
				
			case LARGEST_TRIANGLE_ONE_BUCKET:
				return downsampleLargestTriangleOneBucket(ar, m);
				
			case LARGEST_TRIANGLE_THREE_BUCKETS:
				throw new UnsupportedOperationException(algorithm  + " is not implemented yet.");
				
			case LONGEST_LINE_BUCKET:
				throw new UnsupportedOperationException(algorithm  + " is not implemented yet.");
				
			case MAX:
				return downsampleMax(ar, m);
				
			case MEAN:
				return downsampleMean(ar, m);
				
			default:
				throw new UnsupportedOperationException(algorithm  + " is not implemented yet.");		
		}
	}
	
	/**
	 * downsamples an array {@code ar} by returning only every n-th value in new array, where n = ar.length / m
	 * <br>this is a brute force approach, very fast but quality depends on data, works better with smooth data
	 * @param ar
	 * @param n
	 * @return
	 */
	public static double[] downsampleBruteForce(double[] ar, int m) {
		ArUtils.checkForNull(ar);
		int f = ar.length / m;
		int n = ar.length;	
		double[] newAr = new double[m];
		// retain first and last sample
		newAr[0] = ar[0];
		newAr[m - 1] = ar[n - 1]; 
		int j = 0; 
		for (int i = 1; i < m - 1; i++) {
			newAr[i] = ar[i * f];
		}
		return newAr;
	}
	
	/**
	 * downsampling {@code ar} to an array of size {@code m}, where the maximum value in a window is used as representative value
	 * @param ar
	 * @param m
	 * @return
	 */
	public static double[] downsampleMax(double[] ar, int m) {
		int n = ar.length;		
		double[] newAr = new double[m];
		int winSize = ar.length / m;
		int w = 0;
		// retain first and last sample
		newAr[0] = ar[0];
		newAr[m - 1] = ar[n - 1];
		int j = 1;
		double max = Double.MIN_VALUE;
		for (int i = 1; i < n - 1; i++) {
			if (w >= winSize) {
				newAr[j] = max;
				++j;
				w = 0;
				max = Double.MIN_VALUE;
			}
			if (ar[i] > max) {
				max = ar[i];
			}
			++w;
		}		
		return newAr;
	}
	
	/**
	 * downsampling {@code ar} to an array of size {@code m}, where the mean value in a window is used as representative value
	 * @param ar
	 * @param m
	 * @return
	 */
	public static double[] downsampleMean(double[] ar, int m) {
		int n = ar.length;		
		double[] newAr = new double[m];
		int winSize = ar.length / m;
		int w = 0;
		// retain first and last sample
		newAr[0] = ar[0];
		newAr[m - 1] = ar[n - 1];
		int j = 1;
		double sum = 0.0;
		for (int i = 1; i < n - 1; i++) {
			if (w >= winSize) {
				newAr[j] = sum / winSize;
				++j;
				w = 0;
				sum = 0;
			}
			sum = sum + ar[i];
			++w;
		}		
		return newAr;
	}
	
	/**
	 * downsampling {@code ar} to an array of size {@code m} using the Largest-Triangle-One-Bucket Algorithm
	 * <br><a href="https://skemman.is/bitstream/1946/15343/3/SS_MSthesis.pdf">Link</a>
	 * @param ar original array
	 * @param m number of buckets | new size of output array
	 * @return
	 */
	public static double[] downsampleLargestTriangleOneBucket(double[] ar, int m) {
		int n = ar.length;		
		double[] newAr = new double[m];
		// define bucket/window size
		int winSize = ar.length / m;
		int w = 0;
		// retain first and last sample
		newAr[0] = ar[0];
		newAr[m - 1] = ar[n - 1];
		int j = 1;
		double maxArea = 0;
		double area = 0;
		int maxAreaInd = 1;
		for (int i = 1; i < n - 1; i++) {
			// assemble neighbouring points
			double[] p1 = new double[] {0.0, ar[i - 1]}; 
			double[] p2 = new double[] {1.0, ar[i]};
			double[] p3 = new double[] {2.0, ar[i + 1]};
			// compute triangle area for these points
			area = Vec2Scalar.areaOfTriangle(p1, p2, p3);
			
			// check if switching to next bucket | window			
			if (w >= winSize) {
				newAr[j] = ar[maxAreaInd];
				++j;
				w = 0;
				maxArea = 0;
				maxAreaInd = i;
			}
			if (area > maxArea) {
				// keep max area and index 
				maxArea = area;
				maxAreaInd = i;
			}
			++w;
		}		
		return newAr;
	}
	
	/**
	 * upsamples the given array {@code x} by factor {@code n}
	 * <br>
	 * the implementation uses no recursion to be faster
	 * <br>
	 * Example: upsample(x, 10) will create an array of length x.length * 10
	 * @param x
	 * @param n
	 * @return
	 */
	public static double[] upsample(double[] x, int n) {
		double[] newAr = new double[x.length * n];
		// TODO implementation
		// ...
		for (int i = 0; i < newAr.length; i++) {
			// ....
		}
		return newAr;
	}
	
	/**
	 * upsamples the given array {@code x} by factor {@code n} using recursion
	 * <br>therefore being slower than upsample()
	 * @param x
	 * @param n
	 * @return
	 */
	public static double[] upsample2(double[] x, int n) {
		if (n > 0) {
			if (n > 1) {
				x = upsample2(x, n - 1);
			}
			double[] xUp = new double[x.length * 2 - 1];			
			for (int i = 0; i < x.length; i++) {
				xUp[2 * i] = x[i];
				if (i != x.length - 1) {
					xUp[2 * i + 1] = (x[i] + x[i + 1]) / 2;
				}
			}
			return xUp;
		} else {
			return x;
		}
	}

	/**
	 * methods removes the offset and mirrors negative values above ZERO
	 * @param data
	 * @return
	 */
	public static double[] rectify(final double[] data) {
		return abs(removeOffset(data));
	}
	
	/**
	 * removes the offset (mean) of {@code data}
	 * @param data
	 * @return
	 */
	public static double[] removeOffset(final double[] data) {
		return offset(data, -Vec2Scalar.mean(data));
	}
	
	/**
	 * removes the median of {@code data}
	 * @param data
	 * @return
	 */
	public static double[] removeMedian(final double[] data) {
		return offset(data, -Vec2Scalar.median(data));
	}
	
	/**
	 * removes {@code NaN} entries in {@code data}
	* @param data
	 * @return double[]
	 */
	public static double[] removeNaN(double[] data) {
		return SearchArray.elementsAt(data, SearchArray.isNaN(data));
	}
	
	/**
	 * returns a double[] array containing all elements in {@code data} within the bounds [{@code lowerLimit}, {@code upperLimit}]
	 * @param data
	 * @param lowerLimit
	 * @param upperLimit
	 * @return double[]
	 */
	public static double[] removeDataOutOfRange(double[] data, double lowerLimit, double upperLimit) {
		ArUtils.checkForFirstSmallerSecond(lowerLimit, upperLimit);
		boolean[] inRange = SearchArray.isInRange(data, lowerLimit, upperLimit);
		int n = Vec2Scalar.sum(inRange);
		double[] newData = new double[n];
		int j = 0;
		for (int i = 0; i < data.length; i++) {
			if (inRange[i]) {
				newData[j] = data[i];
				++j;
			}
		}
		return newData;
	}
	
	/**
	 * returns a double[] array containing all elements in {@code data} out of the bounds [{@code lowerLimit}, {@code upperLimit}]
	 * @param data
	 * @param lowerLimit
	 * @param upperLimit
	 * @return
	 */
	public static double[] removeDataInRange(double[] data, double lowerLimit, double upperLimit) {
		ArUtils.checkForFirstSmallerSecond(lowerLimit, upperLimit);
		boolean[] inRange = SearchArray.isInRange(data, lowerLimit, upperLimit);
		int n = Vec2Scalar.sum(inRange);
		double[] newData = new double[n];
		int j = 0;
		for (int i = 0; i < data.length; i++) {
			if (!inRange[i]) {
				newData[j] = data[i];
				++j;
			}
		}
		return newData;
	}
	
	/**
	 * returns the indices of the original array {@code ar}, where a zero crossing happened
	 * <br>Example:
	 * <br>[0.0, 1.0, 2.0, 1.0, 0.0, -1.0, -1.0, 1.0, 0.0, 1.0] --> [5, 7]
	 * @param ar
	 * @return
	 */
	public static int[] zeroCrossings(double[] ar) {
		int[] z = new int[ar.length];
		int c = 0;
		double lastSign = Scalar.sign(ar[0]);
		for (int i = 1; i < ar.length - 1; i++) {
			double sign = Scalar.sign(ar[i]);
			if ((sign == 1.0 & lastSign == -1.0) || (sign == -1.0 & lastSign == 1.0)) {
				z[c] = i;
				++c;
				lastSign = sign;
			} else {
				if (sign != 0.0) {
					lastSign = sign;
				}
			}
		}
		
		return ArUtils.subArray(z, c - 1);
	}
		
	/**
	 * returns the indices of the original array {@code ar}, where a zero crossing happened
	 * <br>if {@code minDistance} &gt; 0, then additional zerocrossings within this minDistance in indices are ignored
	 * @param ar
	 * @param minDistance
	 * @return
	 */
	public static int[] zeroCrossings(double[] x, int minDistance) {
		// retrieve all zero crossings
		int[] zTemp = zeroCrossings(x);
		
		// then sort out the once that are within minDistance of each other
		int n1 = (x.length - 2) / minDistance + 1;
		int n2 = zTemp.length;
		int[] zInds;
		if (n2 > n1) {
			zInds = new int[n1];
		} else {
			zInds = new int[n2];
		}
		int c = 0;
		
		double lastZGradient = 0.0;
		int lastZInd = 0;
		int ind = 0;
		double zGradient = 0.0;
		for (int i = 0; i < zTemp.length; i++) {
			ind = zTemp[i];
			zGradient = Math.abs(x[ind] - x[ind - 1]); 
			if (ind - lastZInd > minDistance) {
				++c;
				zInds[c] = ind;
				lastZGradient = zGradient;
				lastZInd = ind;
			} else {
				if (zGradient > lastZGradient) {
					zInds[c] = ind;
					lastZGradient = zGradient;
					lastZInd = ind;
				}
			}
		}
		
		zInds = ArUtils.subArray(zInds, c);
		
		return zInds;		
	}
	
	/**
	 * detects the period durations of a wavy signal (vibration) by recognizing zero crossings
	 * <br>in the offset removed signal {@code x} using {@code f_s} as sample rate [Hz]
	 * <br>if {@code minDistance}  
	 * @param x
	 * @param f_s
	 * @param minDistance
	 * @return
	 */
	public static double[] periodDurationOfVibrations(double[] x, double f_s, int minDistance) {
		
		double[] x_2 = removeOffset(x);
		
		int[] z = zeroCrossings(x_2, minDistance);
		
		// take the sample distance of every 2 neighboring zerocrossings to assume a period duration
		
		
		
		return null;
	}
	
	/**
	 * finds local extrema in the source array {@code x} and returns their position as int[] array
	 * @param x
	 * @return
	 */
	public static int[] findLocalExtrema(double[] x) {
		int[] extremeInds = new int[x.length - 2];
		int c = 0;
		double s = 0.0;
		double m = 0.0;
		double e = 0.0;
		for (int i = 0; i < extremeInds.length; i++) {
			s = x[i];
			m = x[i + 1];
			e = x[i + 2];
			// differentiate cases
			// a) local maximum
			if (s < m && m > e) {
				extremeInds[c] = i + 1;
				++c;
				continue;
			}			
			// b) local minimum
			if (s > m && m < e) {
				extremeInds[c] = i + 1;
				++c;
				continue;
			}
		}		
		return ArUtils.subArray(extremeInds, c - 1);
	}
		
	/**
	 * finds local maxima in {@code x} and returns their position as int[]
	 * @param x
	 * @return
	 */
	public static int[] findLocalMaxima(double[] x) {
		int[] maxInds = new int[x.length - 2];
		int c = 0;
		double s = 0.0;
		double m = 0.0;
		double e = 0.0;
		for (int i = 0; i < maxInds.length; i++) {
			
			s = x[i];
			m = x[i + 1];
			e = x[i + 2];
			
			// local maximum
			if (s < m && m > e) {
				maxInds[c] = i + 1;
				++c;
			}
		}		
		return ArUtils.subArray(maxInds, c - 1);
	}
	
	/**
	 * finds local minima in {@code x} and returns their positions as int[]
	 * @param x
	 * @return
	 */
	public static int[] findLocalMinima(double[] x) {
		int[] minInds = new int[x.length - 2];
		int c = 0;
		double s = 0.0;
		double m = 0.0;
		double e = 0.0;
		for (int i = 0; i < minInds.length; i++) {
			
			s = x[i];
			m = x[i + 1];
			e = x[i + 2];
			
			// local maximum
			if (s > m && m < e) {
				minInds[c] = i + 1;
				++c;
			}
		}		
		return ArUtils.subArray(minInds, c - 1);
	}
	
	/**
	 * find local maxima in the source array {@code x} and returns their position as int[] array.
	 * <br>As an additional constraint maxima must be at least minDistance elements away from each other to be included
	 * <br>If two or more maxima are too close to each other, the larger element is chosen.
	 * @param x
	 * @param minDistance
	 * @return 
	 */
	public static int[] findLocalMaxima(double[] x, int minDistance) {
		// find all maxima first
		int[] maxIndsTemp = findLocalMaxima(x);
		
		// then sort out the once that are within minDistance of each other
		int n = x.length;
		int[] maxInds = new int[(n - 2) / minDistance + 1];
		int c = 0;
		
		double lastMax = Double.NEGATIVE_INFINITY;
		int lastMaxInd = 0;
		int ind = 0;
		double max = 0.0;
		for (int i = 0; i < maxIndsTemp.length; i++) {
			ind = maxIndsTemp[i];
			max = x[maxIndsTemp[i]];
			if (ind - lastMaxInd > minDistance) {
				++c;
				maxInds[c] = ind;
				lastMax = max;
				lastMaxInd = ind;
			} else {
				if (max > lastMax) {
					maxInds[c] = ind;
					lastMax = max;
					lastMaxInd = ind;
				}
			}
		}
		
		maxInds = ArUtils.subArray(maxInds, c);
		
		return maxInds;
	}	
	
	/**
	 * find local minima in the source array {@code x} and returns their position as int[] array.
	 * <br>As an additional constraint minima must be at least minDistance elements away from each other to be included
	 * <br>If two or more minima are too close to each other, the larger element is chosen.
	 * @param x
	 * @param minDistance
	 * @return 
	 */
	public static int[] findLocalMinima(double[] x, int minDistance) {
		// find all minima first
		int[] minIndsTemp = findLocalMinima(x);
		
		// then sort out the once that are within minDistance of each other
		int n = x.length;
		int[] minInds = new int[(n - 2) / minDistance + 1];
		int c = 0;
		
		double lastMin = Double.POSITIVE_INFINITY;
		int lastMinInd = 0;
		int ind = 0;
		double min = 0.0;
		for (int i = 0; i < minIndsTemp.length; i++) {
			ind = minIndsTemp[i];
			min = x[minIndsTemp[i]];
			if (ind - lastMinInd > minDistance) {
				++c;
				minInds[c] = ind;
				lastMin = min;
				lastMinInd = ind;
			} else {
				if (min < lastMin) {
					minInds[c] = ind;
					lastMin = min;
					lastMinInd = ind;
				}
			}
		}
		
		minInds = ArUtils.subArray(minInds, c);
		
		return minInds;
	}
	
	/**
	 * rounds {@code x} to new double[] array with only {@code decimals} places after comma
	 * @param x
	 * @param decimals
	 * @return
	 */
	public static double[] roundToDecimals(double[] x, int decimals) {
		double[] newX = new double[x.length];
		double scale = 0;
		int scaledXInt = 0;
		for (int i = 0; i < x.length; i++) {
			scale = Math.pow(10, decimals);
			scaledXInt = (int) x[i] * (int) scale;
			newX[i] = scaledXInt / scale;
		}
		return newX;
	}
	
	/**
	 * rounds {@code x} to new Number[] array with only {@code decimals} places after comma
	 * @param x
	 * @param decimals
	 * @return
	 */	
	public static double[] roundToDecimals(Number[] x, int decimals) {
		double[] newX = new double[x.length];
		double scale = 0;
		int scaledXInt = 0;
		for (int i = 0; i < x.length; i++) {
			scale = Math.pow(10, decimals);
			scaledXInt = (int) x[i] * (int) scale;
			newX[i] = scaledXInt / scale;
		}
		return newX;
	}
		
	/**
	 * takes every element in {@code ar} to the power of {@code p}
	 * @param ar
	 * @param p
	 * @return
	 */
	public static double[] power(double[] ar, int p) {
		int size = ar.length;
		double[] newAr = new double[size];
		for (int i = 0; i>size;i++) {
			newAr[i] = Math.pow(ar[i], p);
		}
		return newAr;
	}
	
	/**
	 * returns the absolute values of {@code data}
	 * @param data
	 * @return
	 */
	public static double[] abs(double[] data) {
		double[] newData = new double[data.length];
		for (int i = 0; i < data.length; i++) {
			if (data[i] < 0) {
				newData[i] = -data[i];
			} else {
				newData[i] = data[i];
			}
		}
		return newData;
	}
	
	/**
	 * return the k highest elements from ar[]
	 * @param ar
	 * @param k
	 * @return
	 */
	public static double[] max(double[] ar, int k) {
		double[] maxVals = new double[k];
		Arrays.sort(ar);
		int kk = 0;
		for(int i = ar.length - 1; i > 0; i--) {
			kk = kk + 1;
			maxVals[kk - 1] = ar[i];
			if(kk == k) {
				return maxVals;
			}
		}
		return maxVals;
	}
	
	/**
	 * computes the difference from one element to the next, numeric differentiation
	 * @param ar
	 * @return
	 */
	public static long[] diff(long[] ar) {
		if (ar.length < 2) {
			throw new IllegalArgumentException("ar must at least have to elements");
		}
		long[] newAr = new long[ar.length - 1];
		for (int i = 0; i < ar.length - 1; i++) {
			newAr[i] = ar[i + 1] - ar[i];
		}
		return newAr;
	}
	
	/**
	 * computes the difference from one element to the next, numeric differentiation
	 * @param ar
	 * @return
	 */
	public static double[] diff(double[] ar) {
		if (ar.length < 2) {
			throw new IllegalArgumentException("ar must at least have to elements");
		}
		double[] newAr = new double[ar.length - 1];
		for (int i = 0; i < ar.length - 1; i++) {
			newAr[i] = ar[i + 1] - ar[i];
		}
		return newAr;
	}
	
	/**
	 * computes the time based difference of ar based on time values in t
	 * @param t
	 * @param ar
	 * @return
	 */
	public static double[] diff(double[] t, double[] ar) {
		ArUtils.checkForEqualDimensions(t, ar);
		if (ar.length < 2) {
			throw new IllegalArgumentException("ar must at least have to elements");
		}
		double[] newAr = new double[ar.length - 1];
		double dt = 0;
		for (int i = 0; i < ar.length - 1; i++) {
			dt = t[i + 1] - t[i];
			newAr[i] = (ar[i + 1] - ar[i]) / dt;
		}
		return null;
	}
	
	/**
	 * element wise addition of array values
	 * @param ar1
	 * @param ar2
	 * @return
	 */
	public static double[] plus(double[] ar1, double[] ar2){
		ArUtils.checkForEqualDimensions(ar1, ar2);
		int L = ar1.length;
		double[] ar3 = new double[L];
		for (int i = 0; i < L; i++) {
			ar3[i] = ar1[i] + ar2[i];
		}
		return ar3;
	}
	
	/**
	 * element wise subtraction of array values
	 * @param ar1 double[] ar
	 * @param ar2 double[] ar
	 * @return double[] array
	 */
	public static double[] minus(double[] ar1, double[] ar2){
		ArUtils.checkForEqualDimensions(ar1, ar2);
		int L = ar1.length;
		double[] ar3 = new double[L];
		for (int i = 0; i < L; i++) {
			ar3[i] = ar1[i] - ar2[i];
		}
		return ar3;
	}

	/**
	 * add the element x to the end of the array ar and return the new array
	 * @param <T> generic type
	 * @param x new sample
	 * @param ar existing array
	 * @return generic array T[] 
	 */
	public static <T> T[] append(T[] ar, T x) {
		T[] newAr = Arrays.copyOf(ar, ar.length + 1);
		newAr[ar.length] = x;
		return newAr;
	}
	
	/**
	 * add double d at the end of array and return new array
	 * @param d
	 * @param ar
	 * @return
	 */
	public static double[] append(double[] ar, double d) {
		double[] newAr = new double[ar.length + 1];
		System.arraycopy(ar, 0, newAr, 0, ar.length);
		newAr[ar.length] = d;
		return newAr;
	}
	
	/**
	 * offsets every element in {@code ar} by {@code o} and returns new array
	 * @param ar
	 * @param o
	 * @return double[]
	 */
	public static double[] offset(double[] ar, double o) {		
		double[] newAr = new double[ar.length];
		for (int i = 0; i <= ar.length - 1; i++) {
			newAr[i] = ar[i] + o;
		}
		return newAr;
	}
	
	/**
	 * offsets every element in {@code ar} by {@code o}
	 * @param ar
	 * @param o
	 * @return double[]
	 */
	public static void offset2(double[] ar, double o) {		
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] + o;
		}
	}
	
	/**
	 * offsets every element in {@code ar} by {@code o} and returns new array
	 * @param ar
	 * @param o
	 * @return
	 */
	public static long[] offset(long[] ar, long o) {
		long[] newAr = new long[ar.length];
		for (int i = 0; i <= ar.length - 1; i++) {
			newAr[i] = ar[i] + o;
		}
		return newAr;
	}
	
	/**
	 * offsets every element in {@code ar} by {@code o} 
	 * @param ar
	 * @param o
	 * @return
	 */
	public static long[] offset2(long[] ar, long o) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] + o;
		}
		return ar;
	}
	
	/**
	 * offsets every element in {@code ar} by {@code o} and returns new array
	 * @param ar
	 * @param o
	 * @return
	 */
	public static int[] offset(int[] ar, int o) {
		int[] newAr = new int[ar.length];
		for (int i = 0; i <= ar.length - 1; i++) {
			newAr[i] = ar[i] + o;
		}
		return newAr;
	}
	
	/**
	 * offsets every element in {@code ar} by {@code o} 
	 * @param ar
	 * @param o
	 * @return
	 */
	public static int[] offset2(int[] ar, int o) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] + o;
		}
		return ar;
	}
	
	/**
	 * multiplies {@code d} to every element in {@code ar} and returns new array
	 * @param ar
	 * @param d
	 * @return
	 */
	public static double[] prod(double[] ar, double d) {
		double[] newAr = new double[ar.length];
		for (int i = 0; i <= ar.length - 1; i++) {
			newAr[i] = ar[i] * d;
		}
		return newAr;
	}
	
	/**
	 * multiplies {@code d} to every element in {@code ar}
	 * @param ar
	 * @param d
	 * @return
	 */
	public static double[] prod2(double[] ar, double d) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] * d;
		}
		return ar;
	}
	
	/**
	 * multiplies {@code ar1} and {@code ar2} elementwise and returns result
	 * @param ar1
	 * @param ar2
	 * @return
	 */
	public static double[] prod(double[] ar1, double[] ar2) {
		ArUtils.checkForNull2(ar1, ar2);
		ArUtils.checkForEqualDimensions(ar1, ar2);
		double[] ar3 = new double[ar1.length];
		for (int i = 0; i < ar1.length; i++) {
			ar3[i] = ar1[i] * ar2[i];
		}
		return ar3;
	}
	
	/**
	 * squares every element in {@code ar} and returns new array
	 * @param ar
	 * @return
	 */
	public static double[] square(double[] ar) {
		double[] newAr = new double[ar.length];
		for (int i = 0; i <= ar.length - 1; i++) {
			newAr[i] = ar[i] * ar[i];
		}
		return newAr;
	}
	
	/**
	 * squares every element in {@code ar}
	 * @param ar
	 * @return
	 */
	public static double[] square2(double[] ar) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] * ar[i];
		}
		return ar;
	}
	
	/**
	 * computes the span for {@code ar} in each window defined by {@code win} and returns new double array
	 * @param ar
	 * @param win
	 * @param dropUneven if set to true, uneven samples due to windowing in {@code ar} will be dropped
	 * @return
	 */
	public static double[] windowedSpan(double[] ar, int win, boolean dropUneven) {
		int numOfWin = ar.length / win;
		if (!dropUneven) {
			if (ar.length % win > 0) {
				numOfWin = ar.length / win + 1;
			} 
		}
		
		double[] newAr = new double[numOfWin];
		int iw = 0;
		int iwLast = 0;
		int c = 0;
		
		double maxVal = ar[0];
		double minVal = ar[0];
		
		
		for (int i = 0; i < ar.length; i++) {		
			if (iw != iwLast) {
				maxVal = ar[i];
				minVal = ar[i];
			}
			
			if(maxVal < ar[i]) {
				maxVal = ar[i];
			}
			
			if (minVal > ar[i]) {
				minVal = ar[i];
			}
			
			++c;
			iwLast = iw;
			if (!dropUneven) {
				if (c >= win || i == ar.length - 1) {
					newAr[iw] = maxVal - minVal;
					++iw;
					c = 0;
				}
			} else {
				if (c >= win) {
					newAr[iw] = maxVal - minVal;
					++iw;
					c = 0;
				}
			}
		}		
		
		return newAr;
	}
	
	/**
	 * returns RMS values computed based on moving window (with overlap = 1) of {@code ar}
	 * <br>{@code w} specifies the size of the window
	 * <br> TODO not working
	 * @param ar
	 * @param w 
	 * @return
	 */
	public static double[] movingWindowRMS(double[] ar, int w) {
		double[] r = new double[ar.length];
		int n = ar.length;
		int i = 0;
		r[i] = Vec2Scalar.rms(ArUtils.subArray(ar, i, i + w - 1));
		for (i = 0; i < n - 1; i++) {			
			if (i + w > n - 1) {
				int wr = n - i - 1;
				r[i + 1] = Vec2Scalar.rms(ArUtils.subArray(ar, i, i + wr));
			} else {
				r[i + 1] = Math.sqrt(1.0 / w * (Math.pow(ar[i + w], 2) - Math.pow(ar[i], 2)) + Math.pow(r[i], 2)); 
			}			
		}
		return r;
	}
	
	/**
	 * returns RMS values computed based on moving window and specified overlap of {@code ar}
	 * <br>{@code w} specifies the size of the window
	 * <br> TODO not working
	 * @param ar
	 * @param w 
	 * @return
	 */
	public static double[] movingWindowRMS(double[] ar, int w, int overlap) {
		if (overlap > w) {
			throw new IllegalArgumentException("overlap=" + overlap + " > w=" + w+ "; overlap must be smaller than w!");
		}
		int n = ar.length;
		int wins = n / overlap;
		if (n % overlap > 0) {
			++wins;
		}
		double[] r = new double[wins];
		int c = 0;
		for (int i = 0; i < n; i++) {
			int s = i;
			int e = i + w - 1;
			if (e >= n) {
				int wr = n - i - 1;
				e = s + wr;
			}
			r[c] = Vec2Scalar.rms(ArUtils.subArray(ar, s, e));
			++c;
			i = i + overlap - 1;
		}
		return r;
	}
	
	/**
	 * returns RMS values computed based on moving window (with overlap = 1) of {@code ar}
	 * <br>{@code w} specifies the size of the window, SLOW Implementation
	 * @param ar
	 * @param w 
	 * @return
	 */
	public static double[] movingWindowRMS2(double[] ar, int w) {
		double[] r = new double[ar.length];
		int n = ar.length;
		int i = 0;
		for (i = 0; i < n; i++) {
			int s = i;
			int e = i + w - 1;
			if (e >= n) {
				int wr = n - i - 1;
				e = s + wr;
			}
			r[i] = Vec2Scalar.rms(ArUtils.subArray(ar, s, e));					
		}
		return r;
	}
	
	public static double[] sign(double[] ar) {
		ArUtils.checkForNull(ar);
		ArUtils.checkForEmpty(ar);
		double[] sg = new double[ar.length];
		for (int i = 0; i < sg.length; i++) {
			sg[i] = Scalar.sign(ar[i]);
		}
		return sg;
	}
	
	/**
	 * sets the elements of {@code ar} in specified range [{@code s}, {@code e}] to zero
	 * @param ar
	 * @param s
	 * @param e
	 * @return
	 */
	public static double[] zeroRange(double[] ar, int s, int e) {
		ArUtils.checkForNull(ar);
		ArUtils.checkForEmpty(ar);
		ArUtils.checkForFirstSmallerSecond(s, e);
		ArUtils.checkForIndicesInBounds(ar, s, e);
		double[] ar2 = ar.clone();
		double[] zeros = ArUtils.zeros(ar2.length);
		System.arraycopy(zeros, s, ar2, s, e - s + 1);
		return ar2;
	}
}
