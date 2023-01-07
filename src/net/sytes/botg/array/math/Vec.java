package net.sytes.botg.array.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import net.sytes.botg.array.ConvertArray;
import net.sytes.botg.array.spectrum.WindowFunction;
import net.sytes.botg.array.spectrum.WindowFunction.WindowType;

public class Vec {
	
	private static final long BETTER_OF_AS_STREAM_SIZE = 100_000_000;
	
	public enum DownsamplingAlgorithm {
		BRUTE_FORCE, MAX, MEAN, LARGEST_TRIANGLE_THREE_BUCKETS, LARGEST_TRIANGLE_ONE_BUCKET, LONGEST_LINE_BUCKET
	}
	
	// Suppress default constructor for noninstantiability
	private Vec() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * ----------------------------------------------------------------------------
	 * Vector to Scalar
	 * ----------------------------------------------------------------------------
	 */
	
	/**
	 * retrieves the minimum in {@code x}
	 * @param x
	 * @return
	 */
	public static double min(double[] x) {
		checkForNull(x);
		checkForEmpty(x);
		double minVal = x[0];
		for(double d : x) {
			if(minVal > d) {
				minVal = d;
			} else {
				// do nothing
			}
		}
		return minVal;
	}
	
	/**
	 * returns the span value of {@code x}
	 *<br>Example:
	 *<br>span({1.0, 2.5, 3.4}) --&gt; 3.4 - 1.0 = 2.4 
	 * @param x
	 * @return
	 */
	public static double span(double[] x) {
		double maxVal = x[0];
		double minVal = x[0];
		for(double d : x) {
			if(maxVal > d) {
				// do nothing
			} else {
				maxVal = d;
			}
			if(minVal > d) {
				minVal = d;
			} else {
				// do nothing
			}
		}
		return maxVal - minVal;
	}
	
	/**
	 * retrieves the maximum in {@code x}
	 * @param x
	 * @return
	 */
	public static double max(double[] x) {
		checkForNull(x);
		checkForEmpty(x);
		double maxVal = x[0];
		for(double d : x) {
			if(maxVal > d) {
				// do nothing
			} else {
				maxVal = d;
			}
		}
		return maxVal;
	}	
	
	/**
	 * rms value of {@code x}
	 * @param x
	 * @return
	 */
	public static double rms(double[] x) {
		double rmsSum = 0;
		for (double d : x) {
			rmsSum = rmsSum + d * d;
		}		
		return Math.sqrt(rmsSum / x.length);
	}
	
	/**
	 * rms mean value of {@code x}
	 * @param x
	 * @return
	 */
	public static double rmsMean(double[] x) {
		double rmsSum = 0;
		double sum = 0;
		for (double d : x) {
			rmsSum = rmsSum + d * d;
			sum = sum + d;
		}		
		return Math.sqrt(rmsSum / x.length) + sum / x.length;
	}
	
	/**
	 * sum of {@code x} 
	 * @param x
	 * @return
	 */
	public static double sum(double[] x) {
		double sum = 0;
		for(double d : x) {
			sum = sum + d;
		}
		return sum;
	}
	
	/**
	 * sum of Object Double[] {@code x} 
	 * @param x
	 * @return
	 */
	public static double sum(Double[] x) {
		double sum = 0;
		for(double d : x) {
			sum = sum + d;
		}
		return sum;
	}
	
	/**
	 * sum of {@code x}
	 * @param x
	 * @return
	 */
	public static double sum (int[] x) {
		int sum = 0;
		for(int d : x) {
			sum = sum + d;
		}
		return sum;
	}
	
	/**
	 * sum of {@code x}
	 * @param x
	 * @return
	 */
	public static long sum(long[] x) {
		long sum = 0;
		for(long lo : x) {
			sum = sum + lo;
		}
		return sum;
	}
	
	/**
	 * sum of {@code x} based on Streams
	 * @param ar
	 * @return
	 */
	public static double sum2(double[] x) {
		return Arrays.stream(x).sum();
	}
	
	/**
	 * counts the number of {@code} true elements in {@code x}
	 * @param x
	 * @return int
	 */
	public static int sum(boolean[] x) {
		int sum = 0;
		for (boolean b : x) {
		    sum += b ? 1 : 0;
		}
		return sum;
	}
	
	/**
	 * computes the median
	 * @param x
	 * @return
	 */
	public static double median(double x[]) {
		// array must be cloned before sorting, otherwise the original array is sorted
		double[] x_c = x.clone();
		Arrays.sort(x_c);
		double med = 0;
		if (x_c.length % 2 == 0) {
		    med = (x_c[x_c.length / 2] + x_c[x_c.length / 2 - 1]) / 2;
		} else {
		    med = x_c[x_c.length / 2];
		}
		return med;
	}
	
	/**
	 * computes the mean value of the given {@code x}
	 * @param x
	 * @return double
	 */
	public static double mean(double[] x) {
		return sum(x) / x.length;
	}
	
	/**
	 * computes the mean value of the given {@code x} as Object Primitives Double[]
	 * @param x
	 * @return
	 */
	public static double mean(Double[] x) {
		return sum(x) / x.length;
	}
	
	/**
	 * computes the mean of long Array
	 * @param x
	 * @return
	 */
	public static double mean(long[] x) {
		return (double) sum(x) / x.length;
	}
	
	/**
	 * estimates the variance
	 * @param x
	 * @return
	 */
	public static double variance (double[] x) {
		int size = x.length;
		double mean = mean(x);
		double variance = 0;
		for (double d : x) {
			variance = variance + Math.pow(d - mean, 2);
		}
		return variance / size;
	}
	
	/**
	 * estimates the corrected variance
	 * @param x
	 * @return
	 */
	public static double varianceCorrected(double [] x) {
		int size = x.length;
		double mean = mean(x);
		double variance = 0;
		for (double d : x) {
			variance = variance + Math.pow(d - mean, 2);
		}
		return variance / (size - 1);
	}
	
	/**
	 * estimates the skewness through its moment - biased estimator
	 * @param x
	 * @return
	 */
	public static double skewness(double[] x) {
		int size = x.length;
		double mean = mean(x);
		double var = variance(x);
		double [] x_n = new double [size];
		for (int i =0; i> size; i++) {
			x_n[i] = Math.pow(((x[i] - mean) / var), 3);
		}
		double skewness = mean(x_n);
		return skewness;
	}
	
	/**
	 * estimates the skewness as unbiased estimator
	 * @param x
	 * @return
	 */
	public static double skewnessUnbiased(double[] x) {
		int size = x.length;
		double mean = mean(x);
		double var = variance(x);
		double [] x_n = new double [size];
		for (int i =0; i> size; i++) {
			x_n[i] = Math.pow(((x[i] - mean)/var),3);
		}
		double sum = sum(x_n);
		double skewness = (size/((size-1)*(size-2)))*sum;
		return skewness;
	}
	
	/**
	 * estimates the kurtosis through its moment 
	 * @param x
	 * @return
	 */
	public static double kurtosis(double[] x) {
		int size = x.length;
		double mean = mean(x);
		double var = variance(x);
		double [] x_n = new double [size];
		for (int i =0; i> size; i++) {
			x_n[i] = Math.pow(((x[i] - mean)/var),4);
		}
		double kurtosis = mean(x_n);
		return kurtosis;
	}
	
	/**
	 * computes the crest factor of {@code x}
	 * @param x
	 * @return
	 */
	public static double crest(double[] x) {
		return max(x) / rms(x);
	}
		
	/**
	 * computes the sum product of {@code x1} with {@code x2}
	 * @param x1
	 * @param x2
	 * @return
	 */
	public static double sumprod(double[] x1, double[] x2) {
		checkForEqualDimensions(x1, x2);
		double sumprod = 0;
		for(int i = 0; i < x1.length; i++) {
			sumprod = sumprod + x1[i] * x2[i];
		}
		return sumprod;
	}
	
	/**
	 * computes the scalar product between {@code x1} and {@code x2}
	 * @param x1
	 * @param x2
	 * @return
	 */
	public static double scalarProduct(double[] x1, double[] x2) {
		checkForEqualDimensions(x1, x2);
		double scalar = 0.0;
		for (int i = 0; i < x1.length; i++) {
			scalar += x1[i] * x2[i];
		}
		return scalar;
	}
	
	/**
	 * computes the distance between two points in 1D and greater
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double distance(double[] p1, double[] p2) {
		checkForNull(p1);
		checkForNull(p2);
		checkForEqualDimensions(p1, p2);
		int n = p1.length;
		double sum = 0;
		for (int i = 0; i < n; i++) {
			sum += Math.pow(p1[i] - p2[i], 2);
		}
		return Math.sqrt(sum);
	}
	
	/**
	 * computes the norm (length) of the vector {@code x}
	 * @param x
	 * @return
	 */
	public static double norm(double[] x){
		double squaredSum = 0.0;
		for (int i = 0; i < x.length; i++) {
			squaredSum = squaredSum + x[i] * x[i];
		}
		return Math.sqrt(squaredSum);
	}
	
	/**
	 * computation of area of triangle defined by 3 points [x1, y1], [x2, y2] and [x3, y3]
	 * and Heron' Formula (<a href="https://en.wikipedia.org/wiki/Heron%27s_formula">Link</a>)
	 * <br>
	 * <a href="https://en.wikipedia.org/wiki/File:Triangle_with_notations_2_without_points.svg">Triangle</a>
	 * <br>
	 * @param p1
	 * @param p2
	 * @param p3
	 * @return
	 */
	public static double areaOfTriangle(double[] p1, double[] p2, double[] p3) {
		double a = distance(p1, p2);
		double b = distance(p2, p3);
		double c = distance(p3, p1);
		double s = (a + b + c) / 2;
		double A = Math.sqrt(s * (s - a) * (s - b) * (s - c));
		return A;
	}
	
	/**
	 * computes the mean squared error between two double arrays
	 * @param x1
	 * @param x2
	 * @return
	 */
	public static double mse(double[] x1, double[] x2) {
		checkForEqualDimensions(x1, x2);
		double sum = 0.0;
		for (int i = 0; i < x1.length; i++) {			
			sum = sum + Math.pow(x1[i] - x2[i], 2);
		}
		return sum / x1.length;
	}
	
	/**
	 * computes the mean squared error between two int arrays
	 * @param x2
	 * @param x1
	 * @return 
	 */
	public static double mse(int[] x1, int[] x2) {
		checkForEqualDimensions(x1, x2);
		double sum = 0.0;
		for (int i = 0; i < x1.length; i++) {			
			sum = sum + Math.pow(x1[i] - x2[i], 2);
		}
		return (double) sum / x1.length;
	}
	
	/**
	 * ----------------------------------------------------------------------------
	 * Vector to Vector
	 * ----------------------------------------------------------------------------
	 */
	
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
	 * applies a discrete 1st order low pass filter on the input data {@code x}
	 * <br>specified by resistance {@code R}, capacity {@code C} and sampling period {@code T_s}
	 * <br>returning the filtered signal as double array  
	 * @param x
	 * @param R
	 * @param C
	 * @param T_s
	 * @return
	 */
	public static double[] lowpass(double[] x, double R, double C, double T_s) {
		double[] y = new double[x.length];
		double c = Math.exp(- T_s / R / C);
		y[0] = x[0];
		for (int i = 1; i < y.length; i++) {
			y[i] = (1 - c) * x[i] + c * y[i - 1];
		}
		return y;
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
		checkForNull(x);
		checkForEmpty(y);
		checkForNull(xi);
		checkForEqualDimensions(x, y);
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
	
	public static double[] upperEnvelope(double[] y, int minDistance) {
		int n = y.length;
		double[] x_ue = new double[n];
		double[] x = linspace(n);
		int[] minInds = findLocalMinima(y, minDistance);
		double[] x_s = elementsAt(x, minInds);
		double[] y_s = elementsAt(y, minInds);
		
		return x_ue;
	}
	
	public static double[] lowerEnvelope(double[] x, int minDistance) {
		int n = x.length;
		double[] x_ue = new double[n];
		return x_ue;
	}
	
	/**
	 * returns the y-coordinates of a spline interpolation between the specified supporting points with coordinates {{@code x_support}, {@code y_support}} for the values of {@code x}
	 * <br>the interpolation is done based on cubic natural splines
	 * <br>the following formulas comprise a linear equation system, that can be solved to generate a piece-wise interpolation function for the spline
	 * <br>a) a<sub>i</sub>x<sub>i</sub><sup>3</sup> + b<sub>i</sub>x<sub>i</sub><sup>2</sup> + c<sub>i</sub>x<sub>i</sub> + d<sub>i</sub> = y<sub>i</sub>
	 * <br>b) 3a<sub>i-1</sub>x<sub></sub><sup>2</sup> + 2b<sub>i-1</sub>x<sub>i</sub> + c<sub>i-1</sub> - 3a<sub>i</sub>x<sub>i</sub><sup>2</sup> - 2b<sub>i</sub>x<sub>i</sub> - c<sub>i</sub> = 0
	 * <br>c) 6a<sub>i-1</sub>x<sub>i</sub> + 2b<sub>i-1</sub> - 6a<sub>i</sub>x<sub>i</sub> - 2b<sub>i</sub> = 0, except
	 * <br>d) 6a<sub>1</sub>x<sub>1</sub> + 2b<sub>1</sub> = 0 and 6a<sub>n</sub>x<sub>n</sub> + 2b<sub>n</sub> = 0
	 * @param x_support x-coordinates of supporting points of the spline
	 * @param y_support y-coordinates of supporting points of the spline
	 * @param x x-coordinates the interpolation shall be created for 
	 * @return double[] array of interpolated y-coordinates for {@code x}
	 */
	public static double[] splineInterp(double[] x_support, double[] y_support, double[] x) {
		checkForNull2(x_support, y_support, x);
		checkForEmpty2(x_support, y_support, x);
		checkForEqualDimensions(x_support, y_support);
		
		int order = 3;
		int ns = x_support.length;	// number of supporting points
		int np = ns - 1;	// number of polynoms
		int nr = np * (order + 1);	// number of rows (columns)
		int ni = x.length;	// number of points for interpolation;
		
		double[] y_lin = new double[nr];	// result vector of linear equation system
		double[] x_lin;	// polynom coefficients
		double[][] A = new double[nr][nr];	// linear equation system matrix
		
		double[] y_s = new double[ni];		// y-coordinates of interpolation 
		
		// populate y and A
		// a) Equations expressing that the polynomials traverse the points {x_support, y_support}
		int j = 0;
		int m = 0;
		for (int i = 0; i < ns; i++) {
			if (i == 0) {
				y_lin[j] = y_support[i];
				for (int k = 0; k < order + 1; k++) {
					A[j][k + m * (order + 1)] = Math.pow(x_support[i], order - k);
				}
				++j;
			} else if (i == ns - 1) {
				y_lin[j] = y_support[i];
				for (int k = 0; k < order + 1; k++) {
					A[j][k + m * (order + 1)] = Math.pow(x_support[i], order - k);
				}
				++j;
			} else {
				y_lin[j] = y_support[i];
				for (int k = 0; k < order + 1; k++) {
					A[j][k + m * (order + 1)] = Math.pow(x_support[i], order - k);
				}
				++j;
				++m;
				y_lin[j] = y_support[i];
				for (int k = 0; k < order + 1; k++) {
					A[j][k + m * (order + 1)] = Math.pow(x_support[i], order - k);
				}
				++j;
			}
		}
		
		// b) The first derivative of two adjacent polynomials is equal in the point where they touch
		m = 0;
		for (int i = 1; i < np; i++) {
			y_lin[j] = 0;
			for (int k = 0; k < order + 1; k++) {
				if (order - k <= 0) {
					A[j][k + m * (order + 1)] = 0.0;
					A[j][k + (m + 1) * (order + 1)] = 0.0;
				} else {
					A[j][k + m * (order + 1)] = (order - k) * Math.pow(x_support[i], order - 1 - k);
					A[j][k + (m + 1) * (order + 1)] = -(order - k) * Math.pow(x_support[i], order - 1 - k);
				}
			}
			++m;
			++j;
		}
		
		// c) The second derivative of two adjacent polynomials is equal in the point where they touch
		m = 0;
		for (int i = 1; i < np; i++) {
			y_lin[j] = 0;
			for (int k = 0; k < order + 1; k++) {
				if (order - k - 1 <= 0) {
					A[j][k + m * (order + 1)] = 0.0;
					A[j][k + (m + 1) * (order + 1)] = 0.0;
				} else {
					A[j][k + m * (order + 1)] = (order - k) * (order - k - 1) * Math.pow(x_support[i], order - 2 - k);
					A[j][k + (m + 1) * (order + 1)] = -(order - k) * (order - k - 1) * Math.pow(x_support[i], order - 2 - k);
				}
			}
			++m;
			++j;
		}
		
		// d) The boundary condition (natural spline)
		
		for (int k = 0; k < order + 1; k++) {
			if (order - k - 1 <= 0) {
				A[j][k] = 0.0;
				A[j][k + nr - (order + 1)] = 0.0;
			} else {
				y_lin[j] = 0.0;
				y_lin[j + 1] = 0.0;
				A[j][k] = (order - k) * (order - k - 1) * Math.pow(x_support[0], order - 2 - k);
				A[j + 1][k + nr - (order + 1)] = (order - k) * (order - k - 1) * Math.pow(x_support[ns - 1], order - 2 - k);
			}
		}
		
		// solve the linear equation system
		x_lin = Mat.linSolveGaussian(A, y_lin);
		
		// interpolate the points
		int lastInd = 0;
		for (int i = 0; i < ni; i++) {
			if (x[i] <= x_support[0]) {
				lastInd = 0;
			} else if (x[i] >= x_support[ns - 1]) {
				lastInd = np - 1;
			} else {
				for (int k = lastInd; k < np; k++) {
					if (x[i] >= x_support[k] && x[i] < x_support[k + 1]) {
						lastInd = k;
						break;
					}
				}
			}
			double ps = 0.0;
			for (int k = 0; k < order + 1; k++) {
				ps = ps + x_lin[lastInd * (order + 1) + k] * Math.pow(x[i], order - k);
			}
			y_s[i] = ps;
		}
		
		return y_s;
	}
	
	/**
	 * removes {@code start} samples from the beginning of the {@code ar} and {@code end} samples from the end of {@code ar} 
	 * <br>and returns the result as new double[] array 
	 * @param start
	 * @param end
	 * @return
	 */
	public static double[] trim(double[] ar, int start, int end) {
		checkForGreaterEqualZero(new int[] {start, end});
		int n = ar.length;
		if (start + end > n) {
			throw new IllegalArgumentException("sum of start and end must be smaller than length of ar");
		}
		return sub(ar, start, n - end); 
	}
	
	/**
	 * normalization of vector {@code x} into range [0, 1]
	 * @param x
	 * @return
	 */
	public static double[] normalize(double[] x) {
		double xMax = max(x);
		double xMin = min(x);
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
		double xMax = max(x);
		double xMin = min(x);
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
		double mean = mean(x);
		double sigma = variance(x);
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
	 * @param x
	 * @param n
	 * @return
	 */
	public static double[] downsampleBruteForce(double[] x, int m) {
		checkForNull(x);
		int f = x.length / m;
		int n = x.length;	
		double[] newAr = new double[m];
		// retain first and last sample
		newAr[0] = x[0];
		newAr[m - 1] = x[n - 1]; 
		int j = 0; 
		for (int i = 1; i < m - 1; i++) {
			newAr[i] = x[i * f];
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
			area = Vec.areaOfTriangle(p1, p2, p3);
			
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
		return offset(data, -mean(data));
	}
	
	/**
	 * removes the median of {@code data}
	 * @param data
	 * @return
	 */
	public static double[] removeMedian(final double[] data) {
		return offset(data, -median(data));
	}
	
	/**
	 * removes {@code NaN} entries in {@code data}
	* @param data
	 * @return double[]
	 */
	public static double[] removeNaN(double[] data) {
		return elementsAt(data, isNaN(data));
	}
	
	/**
	 * returns a double[] array containing all elements in {@code data} within the bounds [{@code lowerLimit}, {@code upperLimit}]
	 * @param data
	 * @param lowerLimit
	 * @param upperLimit
	 * @return double[]
	 */
	public static double[] removeDataOutOfRange(double[] data, double lowerLimit, double upperLimit) {
		Scalar.checkForFirstSmallerSecond(lowerLimit, upperLimit);
		boolean[] inRange = isInRange(data, lowerLimit, upperLimit);
		int n = Vec.sum(inRange);
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
		Scalar.checkForFirstSmallerSecond(lowerLimit, upperLimit);
		boolean[] inRange = isInRange(data, lowerLimit, upperLimit);
		int n = sum(inRange);
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
	 * returns the moving average of {@code ar} with a rectangular window of size {@code w}
	 * @param ar
	 * @param window
	 * @return
	 */
	public static double[] movingAverage(double[] ar, int w) {
		int n = ar.length;
		double[] newAr = new double[n];
		for (int i = 0; i < n; i++) {
			
			/* DEBUG
			if (i == 6841) {
				System.out.println(i);
			}
			*/
			
			int e = i + 1;
			int s = e - w;
			if (s < 0) {
				s = 0;
			}
			double sum = 0;
			for (int j = s; j < e; j++) {
				sum = sum + ar[j];
			}
			newAr[i] = sum / (e - s);
		}
		return newAr;
	}
	
	/**
	 * returns the moving average of {@code ar} with a triangular window of size {@code w}
	 * TODO not properly symmetric window
	 * @param ar
	 * @param w
	 * @return
	 */
	public static double[] movingAverageTriangular(double[] ar, int w) {
		int n = ar.length;
		double[] newAr = new double[n];
		double[] triWin = WindowFunction.generate(WindowType.TRIANGULAR, w);
		for (int i = 0; i < n; i++) {					
			int e = i + 1;
			int s = e - w;
			if (s < 0) {
				s = 0;
			}
			double sum = 0;
			if (e - s != w) {
				triWin = WindowFunction.generate(WindowType.TRIANGULAR, e - s);
				if (triWin.length == 1) {
					triWin[0] = 1.0;
				} else if (triWin.length == 2) {
					triWin[0] = 0.5;
					triWin[1] = 0.5;
				}
			}
			int k = 0;
			double wn = 0;
			for (int j = s; j < e; j++) {
				sum = sum + ar[j] * triWin[k];
				wn = wn + triWin[k];
				++k;
			}
			newAr[i] = sum / wn;
		}
		return newAr;		
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
		
		return sub(z, c - 1);
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
		int ind = 0;
		double zGradient = 0.0;
		zInds[0] = zTemp[0];
		int lastZInd = zTemp[0];
		for (int i = 1; i < zTemp.length; i++) {
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
		
		zInds = sub(zInds, c);
		
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
		return sub(extremeInds, c - 1);
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
		return sub(maxInds, c - 1);
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
		return sub(minInds, c - 1);
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
		
		maxInds = sub(maxInds, c);
		
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
		
		minInds = sub(minInds, c);
		
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
		checkForEqualDimensions(t, ar);
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
	 * @param x
	 * @param y
	 * @return
	 */
	public static double[] plus(double[] x, double[] y){
		checkForEqualDimensions(x, y);
		int L = x.length;
		double[] ar3 = new double[L];
		for (int i = 0; i < L; i++) {
			ar3[i] = x[i] + y[i];
		}
		return ar3;
	}
	
	/**
	 * element wise subtraction of array values
	 * @param x double[] array
	 * @param y double[] array
	 * @return double[] array
	 */
	public static double[] minus(double[] x, double[] y){
		checkForEqualDimensions(x, y);
		int L = x.length;
		double[] ar3 = new double[L];
		for (int i = 0; i < L; i++) {
			ar3[i] = x[i] - y[i];
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
	public static double[] product(double[] ar, double d) {
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
	public static double[] product2(double[] ar, double d) {
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
		checkForNull2(ar1, ar2);
		checkForEqualDimensions(ar1, ar2);
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
		r[i] = rms(sub(ar, i, i + w - 1));
		for (i = 0; i < n - 1; i++) {			
			if (i + w > n - 1) {
				int wr = n - i - 1;
				r[i + 1] = rms(sub(ar, i, i + wr));
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
			r[c] = rms(sub(ar, s, e));
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
			r[i] = Vec.rms(sub(ar, s, e));					
		}
		return r;
	}
	
	public static double[] sign(double[] ar) {
		checkForNull(ar);
		checkForEmpty(ar);
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
		checkForNull(ar);
		checkForEmpty(ar);
		Scalar.checkForFirstSmallerSecond(s, e);
		checkForIndicesInBounds(ar, s, e);
		double[] ar2 = ar.clone();
		double[] zeros = zeros(ar2.length);
		System.arraycopy(zeros, s, ar2, s, e - s + 1);
		return ar2;
	}
	
	/**
	 * retrieves the maximum valueand the index it was found in {@code x}
	 * <br>both results are stored in double[], where [0] -&gt; max value and [1] -&gt; index
	 * @param x
	 * @return
	 */
	public static double[] maxI(double[] x) {
		checkForNull(x);
		checkForEmpty(x);
		double maxVal = x[0];
		int maxInd = -1;
		for(int i = 0; i < x.length; i++) {
			if(maxVal > x[i]) {
				// do nothing
			} else {
				maxVal = x[i];
				maxInd = i;
			}
		}
		return new double[] {maxVal, maxInd};
	}	
	
	/**
	 * ----------------------------------------------------------------------------
	 * Vector to Matrix
	 * ----------------------------------------------------------------------------
	 */
	
	/**
	 * reshapes the vector into matrix defined by {@code rows}
	 * <br>the elements are sorted column by column
	 * @param x
	 * @param rows
	 * @return
	 */
	public static double[][] matrix(double[] x, int rows){
		return matrix(x, rows, true);
	}
	
	/**
	 * reshapes the vector into matrix defined by {@code rows}
	 * <br>the elements are sorted column by column if {@code columnByColumn} is set to true
	 * @param x
	 * @param rows
	 * @param columnByColumn
	 * @return
	 */
	public static double[][] matrix(double[] x, int rows, boolean columnByColumn){
		checkForNull(x);
		int n = x.length;
		if (n % rows != 0) {
			throw new IllegalArgumentException("the length of vector x (" + n + ") must be divisble by the number of rows (" + rows + ")");
		}
		int columns = n / rows;
		return matrix(x, rows, columns, columnByColumn);
	}
	
	/**
	 * reshapes the vector into matrix defined by {@code rows} and {@code columns}
	 * <br>the elements are sorted column by column if {@code columnByColumn} is set to true
	 * @param x
	 * @param rows
	 * @param columns
	 * @param columnByColumn
	 * @return
	 */
	public static double[][] matrix(double[] x, int rows, int columns, boolean columnByColumn){
		checkForNull(x);
		checkForEmpty(x);
		int n = x.length;
		if (n != rows * columns) {
			throw new IllegalArgumentException("vector x (" + n + ") must have the same length as rows x columns (" + (rows * columns) + ")");
		}
		double[][] X = new double[rows][columns];
		if (!columnByColumn) {
			for (int i = 0; i < n; i++) {
				int r = i / columns;
				int c = i % columns;
				X[r][c] = x[i];
			}
		} else {
			for (int i = 0; i < n; i++) {
				int c = i / columns;
				int r = i % columns;
				X[r][c] = x[i];
			}
		}
		return X;
	}
	
	/**
	 * separates the given double array {@code ar}
	 * @param ar
	 * @param windowSize
	 * @param overlap
	 * @return
	 */
	public static double[][] overlapWindows(double[] ar, int windowSize, int overlap){
		int n = ar.length;
		if (n < windowSize) {
			throw new IllegalArgumentException("ar.length is smaller than windowSize");
		}
		if (overlap < 0 || overlap > windowSize - 1) {
			throw new IllegalArgumentException("overlap must be in range [0, windowSize - 1]");
		}
		int o;
		int wins;
		int rem;
		if (overlap == 0) {
			wins = n / windowSize;
			rem = n % windowSize;
			if (rem > 0) {
				++wins;
			}
		} else {
			wins = (n - overlap) / (windowSize - overlap);
			rem =  (n - overlap) % (windowSize - overlap);
			if (rem > 0) {
				++wins;
			}
		}
		
		double[][] windows = new double[wins][windowSize];
		int w = 0;
		int j = 0;
		for (int i = 0; i < n; i++) {
			if (j >= windowSize) {
				++w;
				j = 0;
				//i = i - (windowSize - o);
				i = i - overlap;
				if (w == wins - 1) {
					i = n - windowSize;
				}
			}
			windows[w][j] = ar[i];
			++j;
		}
		return windows;
	}
	
	/**
	 * separates the given double array into windows of size windowSize
	 * <br>if dropUneven is true, then the remaining last ar.length % windowSize elements are dropped
	 * <br>dimensions --&gt; double[win][windowSize]
	 * @param ar
	 * @param windowSize
	 * @param dropUneven
	 * @return
	 */
	public static double[][] fixedWindows(double[] ar, int windowSize, boolean dropUneven){
		if (ar == null) {
			return null;
		}
		int win = 0;
		int rem = 0;
		if (dropUneven) {
			win = ar.length / windowSize;
		} else {
			win = ar.length / windowSize;
			rem =  ar.length % windowSize;
			if (rem > 0) {
				++win;
			}
		}
		double[][] wAr = new double[win][];
		for (int i = 0; i < win; i++) {
			double[] dAr = null;
			if (i == win - 1) {
				if (dropUneven) {
					dAr = new double[windowSize];
					System.arraycopy(ar, i * windowSize, dAr, 0, windowSize);
				} else {
					if (rem == 0) {
						dAr = new double[windowSize];
						System.arraycopy(ar, i * windowSize, dAr, 0, windowSize);
					} else {
						dAr = new double[rem];
						System.arraycopy(ar, i * windowSize, dAr, 0, rem);
					}
				}
			} else {
				dAr = new double[windowSize];
				System.arraycopy(ar, i * windowSize, dAr, 0, windowSize);
			}
			wAr[i] = dAr;
		}
		return wAr;
	}

	/**
	 * returns the array {@code data} split into windows of size {@code window}
	 * if dropUneven is set to {@code false}, then the remaining elements are filled with NaN
	 * @param data
	 * @param window
	 * @param dropUneven
	 * @return double[][]
	 */
	public static double[][] separateDataIntoWindows(double[] data, int window, boolean dropUneven){
		double[][] dataWindows = null;
		if (dropUneven) {
			int w = data.length / window;
			dataWindows = new double[w][window];
			int s = 0;
			int e = window;
			for (int i = 0; i < w; i++) {
				s = window * i;
				System.arraycopy(data, s, dataWindows[i], 0, e);
			}
		} else {
			int w = data.length / window;
			int rem = data.length % window;
			if (rem > 0) {
				 ++w;
			}
			dataWindows = new double[w][window];
			int s = 0;
			int e = window -1;
			for (int i = 0; i < w; i++) {
				if (i == w - 1 && rem > 0) {
					double[] nanAr = nan(rem);
					s = window * i;
					System.arraycopy(data, s, dataWindows[i], 0, rem);
					System.arraycopy(nanAr, 0, dataWindows[i], rem, e - rem + 1);
				} else {
					s = window * i;
					System.arraycopy(data, s, dataWindows[i], 0, e);
				}
			}			
		}
		return dataWindows;
	}
	
	/**
	 * ----------------------------------------------------------------------------
	 * Vector Permutation Methods
	 * ----------------------------------------------------------------------------
	 */
	
	public static double[] bubbleSort(double[] ar) {
		for (int a = ar.length - 1; a > 0; a--) {
			for (int b = 0; b < a; b++) {
				if (ar[b] > ar[b + 1]) {
					double temp = ar[b];
					ar[b] = ar[b + 1];
					ar[b + 1] = temp;
				}
			}
		}
		
		return ar;
	}
	
	/**
	 * subject to EPL 2.0 in /lic/LICENSE_VOGELLA_QUICKSORT.txt
	 * @see <a href="https://www.vogella.com/tutorials/JavaAlgorithmsQuicksort/article.html">Link</a>
	 * @param ar
	 */
	public static void quicksort(double[] ar) {		
        double[] numbers;
        int number;
		// check for empty or null array
        if (ar == null || ar.length == 0){
            return;
        }
        numbers = ar;
        number = ar.length;
        quicksortRec(0, number - 1, numbers);
    }
	
	/**
	 * returns a sorted double array based on the elements in {@code ar}
	 * subject to EPL 2.0 in /lic/LICENSE_VOGELLA_QUICKSORT.txt
	 * @see <a href="https://www.vogella.com/tutorials/JavaAlgorithmsQuicksort/article.html">Link</a>
	 * @param ar
	 * @return
	 */
	public static int[] quicksort2(double[] ar) {
        double[] numbers;
        int[] inds = new int[ar.length];
        int number;
		// check for empty or null array
        if (ar == null || ar.length == 0){
            return null;
        }
        numbers = ar;
        number = ar.length;
        for(int i = 0; i < inds.length; i++) {
        	inds[i] = i;
        }
        quicksortRec2(0, number - 1, numbers, inds);
        return inds;
    }

    private static void quicksortRec2(int low, int high, double[] numbers, int[] inds) {
    	int i = low, j = high;
        // Get the pivot element from the middle of the list
        double pivot = numbers[low + (high-low)/2];

        // Divide into two lists
        while (i <= j) {
            // If the current value from the left list is smaller than the pivot
            // element then get the next element from the left list
            while (numbers[i] < pivot) {
                i++;
            }
            // If the current value from the right list is larger than the pivot
            // element then get the next element from the right list
            while (numbers[j] > pivot) {
                j--;
            }

            // If we have found a value in the left list which is larger than
            // the pivot element and if we have found a value in the right list
            // which is smaller than the pivot element then we exchange the
            // values.
            // As we are done we can increase i and j
            if (i <= j) {
            	exchangeQuickSortElements2(i, j, numbers, inds);
                i++;
                j--;
            }
        }
        // Recursion
        if (low < j) {
            quicksortRec2(low, j, numbers, inds);
        }
        if (i < high) {
        	quicksortRec2(i, high, numbers, inds);
        }
	}

	private static void exchangeQuickSortElements2(int i, int j, double[] numbers, int[] inds) {
		// TODO Auto-generated method stub
		double temp = numbers[i];
        numbers[i] = numbers[j];
        numbers[j] = temp;
        int tempInd = inds[i];
        inds[i] = inds[j];
        inds[j] = tempInd;
	}

	private static void quicksortRec(int low, int high, double[] numbers) {
        int i = low, j = high;
        // Get the pivot element from the middle of the list
        double pivot = numbers[low + (high-low)/2];

        // Divide into two lists
        while (i <= j) {
            // If the current value from the left list is smaller than the pivot
            // element then get the next element from the left list
            while (numbers[i] < pivot) {
                i++;
            }
            // If the current value from the right list is larger than the pivot
            // element then get the next element from the right list
            while (numbers[j] > pivot) {
                j--;
            }

            // If we have found a value in the left list which is larger than
            // the pivot element and if we have found a value in the right list
            // which is smaller than the pivot element then we exchange the
            // values.
            // As we are done we can increase i and j
            if (i <= j) {
            	exchangeQuickSortElements(i, j, numbers);
                i++;
                j--;
            }
        }
        // Recursion
        if (low < j) {
            quicksortRec(low, j, numbers);
        }
        if (i < high) {
        	quicksortRec(i, high, numbers);
        }
    }

    private static void exchangeQuickSortElements(int i, int j, double[] numbers) {
        double temp = numbers[i];
        numbers[i] = numbers[j];
        numbers[j] = temp;
    }	
    
	public static void flip(byte[] array) {
		//long start = System.nanoTime();
		int i = 0;
		byte b;
		int n = array.length;
		for (i = 0; i < n / 2; i++) {
			b = array[i];
			array[i] = array[n - i - 1];
			array[n - i - 1] = b;
		}
		//System.out.println(System.nanoTime() - start + " ns elapsed");
	}
	
	// DEBUGGING NOT WORKING DUE TO LIST IMPLEMENTATION
	public static void flip2(byte[] array) {
		//long start = System.nanoTime();
		List<byte[]> bL = Arrays.asList(array);
		Collections.reverse(bL);
		//System.out.println(System.nanoTime() - start + " ns elapsed");
	}

	/**
	 * 
	 * @param array
	 * @return
	 */
	public static byte[] flip3(byte[] array) {
		//long start = System.nanoTime();
		int i, j, n;
		n = array.length;
		j = n;
		byte[] newArray = new byte[j];
		for (i = 0; i < n; i++) {
			newArray[j - 1] = array[i];
			j = j - 1;
		}
		return newArray;
		//System.out.println(System.nanoTime() - start + " ns elapsed");
	}
	
	/**
	 * switches values at odd indices to even indices and vice versa
	 * <br>faster than oddToEven
	 * @param array
	 * @return
	 */
	public static byte[] oddToEven(byte[] array) {
		//long start = System.nanoTime();
		byte[] newArray = new byte[array.length];
		if (array.length < 2) {
			throw new RuntimeException("array must contain at least 2 elements");
		}
		for (int i = 1; i < array.length; i += 2) {
			newArray[i - 1] = array[i];
			newArray[i] = array[i - 1];
		}
		return newArray;
		//System.out.println(System.nanoTime() - start + " ns elapsed");
	}
	
	/**
	 * switches values at odd indices to even indices and vice versa 
	 * @param array
	 */
	public static void oddToEven2(byte[] array) {
		if (array.length < 2) {
			throw new RuntimeException("array must contain at least 2 elements");
		}
		byte b1 = 0;
		byte b2 = 0;
		for (int i = 1; i < array.length; i += 2) {
			b1 = array[i];
			b2 = array[i - 1];
			array[i - 1] = b1;
			array[i] = b2;
		}
	}
	
	/**
	 * switches values at odd indices to even indices and vice versa
	 * @param array
	 * @return
	 */
	public static double[] oddToEven(double[] array) {
		//long start = System.nanoTime();
		double[] newArray = new double[array.length];
		if (array.length < 2) {
			throw new RuntimeException("array must contain at least 2 elements");
		}
		for (int i = 1; i < array.length; i += 2) {
			newArray[i - 1] = array[i];
			newArray[i] = array[i - 1];
		}
		return newArray;
		//System.out.println(System.nanoTime() - start + " ns elapsed");
	}
	
	/**
	 * returns only array values at even indices of d 
	 * @param d
	 * @return
	 */
	public static double[] even(double[] d) {
		int n = d.length / 2 + 1;
		double[] e = new double[n];
		for (int k = 0; k < n; k++) {
            e[k] = d[2*k];
        }
		return e;
	}
	
	/**returns only array values at odd indices of d
	 * 
	 * @param d
	 * @return
	 */
	public static double[] odd(double[] d) {
		int n = d.length / 2;
		double[] o = new double[n];
		for (int k = 0; k < n; k++) {
            o[k] = d[2 * k + 1];
        }
		return o;
	}
	
	/**
	 * flips the passed double array
	 * <br>last element is new first element, and so on
	 * <br>Example:
	 * <br>[1, 2, 3, 4, 5] -&gt; [5, 4, 3, 2, 1]
	 * <br>
	 * <br>this code is based on org.apache.commons.commons-lang3.ArrayUtils.reverse()
	 * @param ar
	 */
	public static void flip(double[] ar) {
		double tmp;
		int i = 0;
		int j = ar.length;
        while (j > i) {
            tmp = ar[j];
            ar[j] = ar[i];
            ar[i] = tmp;
            j--;
            i++;
        }
	}
	
	/**
	 * flips the passed intArray
	 * <br>last element is new first element, and so on
	 * <br>Example:
	 * <br>[1, 2, 3, 4, 5] -&gt; [5, 4, 3, 2, 1]
	 * <br>
	 * <br>this code is based on org.apache.commons.commons-lang3.ArrayUtils.reverse()
	 * @param ar
	 */
	public static void flip(int[] ar) {
		int tmp;
		int i = 0;
		int j = ar.length;
        while (j > i) {
            tmp = ar[j];
            ar[j] = ar[i];
            ar[i] = tmp;
            j--;
            i++;
        }
	}
	
	/**
	 * flips the passed longArray
	 * <br>last element is new first element, and so on
	 * <br>Example:
	 * <br>[1, 2, 3, 4, 5] -&gt; [5, 4, 3, 2, 1]
	 * <br>
	 * <br>this code is based on org.apache.commons.commons-lang3.ArrayUtils.reverse()
	 * @param ar
	 */
	public static void flip(long[] ar) {
		long tmp;
		int i = 0;
		int j = ar.length;
        while (j > i) {
            tmp = ar[j];
            ar[j] = ar[i];
            ar[i] = tmp;
            j--;
            i++;
        }
	}
	
	/**
	 * flips the passed shortArray
	 * <br>last element is new first element, and so on
	 * <br>Example:
	 * <br>[1, 2, 3, 4, 5] -&gt; [5, 4, 3, 2, 1]
	 * <br>
	 * <br>this code is based on org.apache.commons.commons-lang3.ArrayUtils.reverse()
	 * @param ar
	 */
	public static void flip(short[] ar) {
		short tmp;
		int i = 0;
		int j = ar.length;
        while (j > i) {
            tmp = ar[j];
            ar[j] = ar[i];
            ar[i] = tmp;
            j--;
            i++;
        }
	}
	
	/**
	 * ----------------------------------------------------------------------------
	 * Vector Search Methods
	 * ----------------------------------------------------------------------------
	 */
	
	/**
	 * method returns true/false whether any of the array elements ar are contained in String s
	 * @param s
	 * @param ar
	 * @return 
	 */
	public static boolean contains(double[] x, double d) {
		for (double xi : x) {
			if (xi == d) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * method returns true/false whether any of the array elements ar are contained in String s
	 * @param s
	 * @param ar
	 * @return 
	 */
	public static boolean contains(int[] x, int i) {
		for (int xi : x) {
			if (xi == i) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * method returns true/false whether any of the array elements {@code x} contain String {@code s}
	 * @param s
	 * @param ar
	 * @return 
	 */
	public static boolean contains(String[] x, String s) {
		for (String xi : x) {
			if (xi.contains(s)) {
				return true;
			}
		}
		return false;
	}
	
    public static int countOccurences(String[] ar, String s, boolean caseSensitive) {
    	int count = -1;
        for (int i = 0; i < ar.length; i++) {
        	if (!caseSensitive) {
	        	if (s.toLowerCase().equals(ar[i].toLowerCase())){
	        		count++;
	            }
        	} else {
        		if (s.equals(ar[i])){
	        		count++;
	            }
        	}
        }
        return count;
    }   
    
    public static boolean isIntInArray(int[] ar, int i) {
    	for (int j = 0; j < ar.length; j++) {
    		if (ar[j] == i) {
    			return true;
    		}
    	}
    	return false;
    }
	
	/**
	 * returns a sub array starting at index s and ending at index e from ar
	 * @param ar
	 * @param s
	 * @param e
	 * @return
	 */
	public static double[] sub(final double[] ar, int s, int e) {
		checkForIndicesInBounds(ar, s, e);
		double[] ar2 = new double[e - s + 1];
		System.arraycopy(ar, s, ar2, 0, ar2.length);
		return ar2;
	}
	
	/**
	 * returns a sub array starting at 0 and ending at index e from ar
	 * @param ar
	 * @param e
	 * @return
	 */
	public static double[] sub(double[] ar, int e) {
		return sub(ar, 0, e);
	}
	
	/**
	 * returns a sub array starting at 0 and ending at index e from ar
	 * @param ar
	 * @param e
	 * @return
	 */
	public static int[] sub(int[] ar, int e) {
		return sub(ar, 0, e);
	}
	
	/**
	 * returns a sub array starting at index s and ending at index e from ar
	 * @param ar
	 * @param s
	 * @param e
	 * @return
	 */
	public static int[] sub(final int[] ar, int s, int e) {
		checkForIndicesInBounds(ar, s, e);
		int[] ar2 = new int[e - s + 1];
		System.arraycopy(ar, s, ar2, 0, ar2.length);
		return ar2;
	}
	
	/**
	 * returns a sub array starting at index s and ending at index e from ar
	 * @param ar
	 * @param s
	 * @param e
	 * @return
	 */
	public static String[] sub(final String[] ar, int s, int e) {
		checkForIndicesInBounds(ar, s, e);
		String[] ar2 = new String[e - s + 1];
		System.arraycopy(ar, s, ar2, 0, ar2.length);
		return ar2;
	}
	
	/**
	 * return the indices of the k highest elements from ar[]
	 * @param ar
	 * @param k
	 * @return
	 */
	public static int[] maxInd(double[] ar, int k) {
		int[] sortInds = quicksort2(ar);
		int[] maxInds = new int[k];
		int kk = 0;
		for(int i = sortInds.length - 1; i > 0; i--) {
			kk = kk + 1;
			maxInds[kk - 1] = sortInds[i];
			if(kk == k - 1) {
				return maxInds;
			}
		}
		return null;
	}
	
	/**
	 * searches the first index in {@code x}, where the element of {@code x} is greater or equal than {@code d}
	 * <br>if equalOrGreater is set to true, the comparison is >= instead of >
	 * <br>if no element is found -1 is returned
	 * @param x
	 * @param d
	 * @param equalOrGreater
	 * @return
	 */
	public static int findFirstValueGreaterThan(double[] x, double d, boolean equalOrGreater) {
		for (int i = 0; i < x.length; i++) {
			if (equalOrGreater) {
				if (x[i] >= d) {
					return i;
				}
			} else {
				if (x[i] > d) {
					return i;
				}
			}
		}
		return -1;
	}
	
	/**
	 * searches the first index in {@code x}, where the element of {@code x} is smaller or equal than {@code d}
	 * <br>if equalOrSmaller is set to true, the comparison is <= instead of <
	 * <br>if no element is found -1 is returned
	 * @param x
	 * @param d
	 * @param equalOrSmaller
	 * @return
	 */
	public static int findFirstValueSmallerThan(double[] x, double d, boolean equalOrSmaller) {
		for (int i = 0; i < x.length; i++) {
			if (equalOrSmaller) {
				if (x[i] <= d) {
					return i;
				}
			} else {
				if (x[i] < d) {
					return i;
				}
			}
		}
		return -1;
	}

	/**
	 * searches the first index in {@code x}, where the element of {@code x} is equal to {@code d}
	 * <br>if no element is found -1 is returned
	 * @param x
	 * @param d
	 * @return
	 */
	public static int findValue(double[] x, double d) {
		for (int i = 0; i < x.length; i++) {
			if (x[i] == d) {
				return i;
			}
		}
		return -1;
	}
	
	/**
	 * returns an array with unique elements based on values in {@code x}
	 * @param x
	 * @return
	 */
	public static double[] unique(final double[] x) {
		checkForNull(x);
		checkForEmpty(x);
		int n = x.length;
		//HashMap<Double, Double> xu = new LinkedHashMap<Double, Double>();
		HashMap<Double, Double> xu = new HashMap<Double, Double>();
		for (int i = 0; i < n; i++) {
			xu.put(x[i], x[i]);
		}
		return ConvertArray.unwrap(xu.values().toArray(new Double[xu.size()]));
	}
	
	/**
	 * returns a boolean[] array containing true for each element in {@code data} that is within the bounds [{@code lowerLimit}, {@code upperLimit}]
	 * @param data
	 * @param lowerLimit
	 * @param upperLimit
	 * @return
	 */
	public static boolean[] isInRange(double[] data, double lowerLimit, double upperLimit) {
		boolean[] inds = new boolean[data.length];
		for (int i = 0; i < data.length; i++) {
			if (data[i] >= lowerLimit && data[i] <= upperLimit) {
				inds[i] = true;
			}
		}
		return inds;
	}
	
	/**
	 * returns the elements specified with {@code inds} in {@code data} into new array
	 * @param inds
	 * @param data
	 * @return double[]
	 */
	public static double[] elementsAt(double[] data, boolean[] inds) {
		if (data.length != inds.length) {
			throw new IllegalArgumentException("length of data and inds array must have same length");
		}
		int n = Vec.sum(inds);
		double[] newData = new double[n];
		int j = 0;
		for (int i = 0; i < data.length; i++) {
			if (inds[i]) {
				newData[j] = data[i];
				++j;
			}
		}
		return newData;
	}
	
	/**
	 * returns the elements at indices {@code inds} in {@code ar} into new array
	 * @param inds
	 * @param data
	 * @return double[]
	 */
	public static double[] elementsAt(double[] ar, int[] inds) {
		int n = inds.length;
		double[] ar2 = new double[n];
		for (int i = 0; i < n; i++) {
			ar2[i] = ar[inds[i]];
		}
		return ar2;
	}
	
	/**
	 * returns an boolean[] array that is true for every element in {@code data} with {@code NaN}
	 * @param data
	 * @return boolean[]
	 */
	public static boolean[] isNaN(double[] data) {
		boolean[] inds = new boolean[data.length];
		for (int i = 0; i < data.length; i++) {
			// returns true only if the element is NaN
			if (data[i] != data[i]) {
				inds[i] = true;
			}
		}
		return inds;
	}
	
	/**
	 * returns the index of {@code str} found in {@code ar}
	 * <br>if {@code ar} does not contain {@code str}, then -1 is returned 
	 * @param ar
	 * @param str
	 * @return
	 */
	public static int indexAt(final String[] ar, final String str) {
		int i = -1;
		for (int s = 0; s < ar.length; s++) {
			if (ar[s].contentEquals(str)) {
				i = s;
				break;
			}
		}
		return i;
	}
	
	/**
	 * returns the index of {@code d} found in {@code x}
	 * <br>if {@code x} does not contain {@code d}, then -1 is returned 
	 * @param x
	 * @param d
	 * @return
	 */
	public static int indexAt(final double[] x, final double d) {
		int i = -1;
		for (int j = 0; j < x.length; j++) {
			if (x[j] == d) {
				i = j;
				break;
			}
		}
		return i;
	}
	
	public static <T> T[] removeEmpty(final T[] ar) {
		List<T> list = new ArrayList<T>();
		for (int a = 0; a < ar.length; a++) {
			if (ar[a] != null) {
				list.add(ar[a]);
			}
		}
		return list.toArray(Arrays.copyOf(ar, list.size()));
	}	

	public static <T> T[] removeWith(T[] ar, T removeContent) {
		List<T> list = new ArrayList<T>();
		for (int a = 0; a < ar.length; a++) {
			if (!ar[a].equals(removeContent)) {
				list.add(ar[a]);
			}
		}
		return list.toArray(Arrays.copyOf(ar, list.size()));
	}
		
	/**
	 * ----------------------------------------------------------------------------
	 * Vector Generation Methods
	 * ----------------------------------------------------------------------------
	 */
		
	public static int[] randInt(int n, int min, int max) {
		if (min >= max) {
			throw new IllegalArgumentException("max must be greater than min");
		}
		int[] rInts = new int[n];		
		Random r = new Random();
		for (int i = 0; i < n; i++) {
			rInts[i] = r.nextInt((max - min) + 1) + min;
		}
		return rInts;
	}
	
	public static int[] randInt(int size) {
		int[] data;
		if (size < BETTER_OF_AS_STREAM_SIZE) {
			data = new int[size];
			Random r = new Random();
			for (int i = 0; i < size; i++) {
				data[i] = r.nextInt();
			}
		} else {
			data = new Random().ints(size).toArray();
		}
		return data;
	}	
	
	/**
	 * returns an array of size n with 0's
	 * @param n
	 * @return
	 */
	public static double[] zeros(int n) {
		return new double[n];
	}
	
	/**
	 * returns an float array of size n with 0's
	 * @param n
	 * @return
	 */
	public static float[] zerosF(int n) {
		return new float[n];
	}

	/**
	 * returns an array of size {@code n} with 1's
	 * @param n
	 * @return
	 */
	public static double[] ones(int n) {
		double[] ar = new double[n];
		for (int i = 0; i < n; i++) {
			ar[i] = 1;
		}
		return ar;
	}
	
	/**
	 * returns an float array of size {@code n} with 1's
	 * @param n
	 * @return
	 */
	public static float[] onesF(int n) {
		float[] x = new float[n];
		for (int i = 0; i < n; i++) {
			x[i] = 1;
		}
		return x;
	}
	
	/**
	 * returns an int array of size {@code n} with 1's
	 * @param n
	 * @return
	 */
	public static int[] onesI(int n) {
		int[] x = new int[n];
		for (int i = 0; i < n; i++) {
			x[i] = 1;
		}
		return x;
	}
	
	/**
	 * creates a new random double array with values in [0, 1]
	 * @param n
	 * @return
	 */
	public static double[] rand(int n) {
		Random r = new Random();
		double[] ar = new double[n];
		for (int i = 0; i < n; i++) {
			ar[i] = r.nextDouble();
		}
		return ar;
	}
	
	/**
	 * creates a new random double array with values in [{@code start},{@code end}]
	 * @param n
	 * @param start
	 * @param end
	 * @return
	 */
	public static double[] rand(int n, double start, double end) {
		Scalar.checkForFirstSmallerSecond(start, end);
		Random r = new Random();
		double[] ar = new double[n];
		for (int i = 0; i < n; i++) {
			ar[i] = (end - start) * r.nextDouble() + start;
		}
		return ar;
	}
	
	/**
	 * returns a double array starting from {@code start} to {@code end} with {@code step}
	 * @param start
	 * @param end
	 * @param step
	 * @return
	 */
	public static double[] linspace(double start, double end, double step) {
		int size = (int) ((end - start ) / step);
		double[] ar = new double[size];
		for (int i = 0; i < size; i++) {
			ar[i] = start + step * i;
		}
		return ar;
	}
	
	/**
	 * returns a int array starting from {@code start} to {@code end} with {@code step}
	 * @param start
	 * @param end
	 * @param step
	 * @return
	 */
	public static int[] linspace(int start, int end, int step) {
		int size = (int) ((end - start ) / step);
		int[] ar = new int[size];
		for (int i = 0; i < size; i++) {
			ar[i] = start + step * i;
		}
		return ar;
	}
	
	/**
	 * returns a double array starting from {@code start} to {@code end} with equally spaced steps, so that the {@code size} is met
	 * @param start
	 * @param end
	 * @param size
	 * @return
	 */
	public static double[] linspace(double start, double end, int size) {
		double[] ar = new double[size];
		double step = (end - start) / (size - 1);
		for (int i = 0; i < size; i++) {
			ar[i] = start + step * i;
		}
		return ar;
	}
	
	/**
	 * returns a int array starting from {@code start} to {@code end} with step size 1
	 * @param start
	 * @param end
	 * @return
	 */
	public static int[] linspace(int start, int end) {
		int size = (int) ((end - start ) / 1);
		int[] ar = new int[size];
		for (int i = 0; i < size; i++) {
			ar[i] = start + 1 * i;
		}
		return ar;
	}
	
	/**
	 * returns a double array starting from {@code start} with {@code size} steps of {@code step}
	 * @param start
	 * @param size
	 * @param step
	 * @return
	 */
	public static double[] linspace(double start, int size, double step) {
		double[] ar = new double[size];
		ar[0] = start;
		for (int i = 1; i < size; i++) {
			ar[i] = ar[i - 1] + step;
		}
		return ar;
	}
	
	/**
	 * returns a long array starting from {@code start} with size {@code size} and steps of {@code step}
	 * @param start
	 * @param size
	 * @param step
	 * @return
	 */
	public static long[] linspace(long start, int size, long step) {
		long[] ar = new long[size];
		ar[0] = start;
		for (int i = 1; i < size; i++) {
			ar[i] = ar[i - 1] + step;
		}
		return ar;
	}
	
	/**
	 * returns a double array starting from 0 and incrementing by 1 for {@code size} times
	 * @param size
	 * @return
	 */
	public static double[] linspace(int size) {
		return linspace(0, size, 1.0);
	}
	
	/**
	 * Create vector of {@code n} logarithmically spaced values between {@code start} and {@code end} 
	 * @param start
	 * @param end
	 * @param n
	 * @return
	 */
	public static double[] logspace(double start, double end, int n) {
		double[] ar = new double[n];
		double step = (end - start) / (n - 1);
		for (int i = 0; i < n; i++) {
			ar[i] = Math.pow(start + step * i, 10);
		}
		return ar;
	}

	/**
	 * returns a vector with NaN values of length {@code n}
	 * @param size
	 * @return
	 */
	public static double[] nan(int n) {
		double[] data = new double[n];
		for (int i = 0; i < n; i++) {
			data[i] = Double.NaN;
		}
		return data;
	}
	
	/**
	 * ----------------------------------------------------------------------------
	 * Vector Utility Methods
	 * ----------------------------------------------------------------------------
	 */
	
	/**
	 * TODO testing
	 * appends all specified {@code arrays} to a new array consisting of all elements
	 * @param arrays
	 * @return
	 */
	public static double[] append(final double[] ... arrays) {
	    int size = 0;
	    for (double[] a: arrays) {
	        size += a.length;
	        double[] res = new double[size];
	        int destPos = 0;
	        for ( int i = 0; i < arrays.length; i++ ) {
	            if ( i > 0 ) {
	            	destPos += arrays[i-1].length;
	            }
	            int length = arrays[i].length;
	            System.arraycopy(arrays[i], 0, res, destPos, length);
	        }
	        return res;
	    }
		return null;
	}
	
	/**
	 * copies the specified array {@code x}
	 * @param x
	 * @return
	 */
	public static double[] copy(double[] x) {
		return x.clone();
	}
	
	/**
	 * returns true if set specified object {@code o} is an array
	 * @param o
	 * @return
	 */
	public static boolean isArray(Object o) {
		if (o != null && o.getClass().isArray()) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * prints a 1D Object[] array
	 * @param ar
	 */
	public static void print(final Object[] ar) {
		System.out.println(Arrays.toString(ar));
	}
	
	/**
	 * prints a 1D double[] array
	 * @param ar
	 */
	public static void print(final double[] ar) {
		System.out.println(Arrays.toString(ar));
	}
	
	/**
	 * prints a 1D int[] array
	 * @param ar
	 */
	public static void print(final int[] ar) {
		System.out.println(Arrays.toString(ar));
	}
	
	public static void checkForEmpty(double[] ar) {
		if (ar.length == 0) {
			throw new IllegalArgumentException("ar must not be empty");
		}
	}
	
	public static void checkForEmpty2(double[] ... x) {
		for (int i = 0; i < x.length; i++) {
			if (x[i].length == 0) {
				throw new IllegalArgumentException("input arrays must not be empty");
			}
		}
	}	
		
	/**
	 * checks if the given {@code inds} are greater than 0
	 * @param inds
	 */
	public static void checkForGreaterZero(final int ... inds) {
		for (int i : inds) {
			if (i < 0) {
				throw new IllegalArgumentException("Indices must be greater than 0.");
			}
		}
	}
	
	/**
	 * checks if the given {@code inds} are greater than or equal to 0
	 * <br>(the input array is final)
	 * @param inds
	 */
	public static void checkForGreaterEqualZero(final int[] inds) {
		for (int i : inds) {
			if (i <= 0) {
				throw new IllegalArgumentException("Indices must be greater than or equal to 0.");
			}
		}
	}
	
	public static void checkForGreaterEqualZero2(int[] inds) {
		for (int i : inds) {
			if (i <= 0) {
				throw new IllegalArgumentException("Indices must be greater than or equal to 0.");
			}
		}
	}
	
	/**
	 * checks if the start and end indices are within array bounds
	 * @param ar
	 * @param s
	 * @param e
	 */
	public static void checkForIndicesInBounds(double[] ar, int s, int e) {
		if (s < 0 || e >= ar.length) {
			throw new IndexOutOfBoundsException("Index out of Bounds, s=" + s + " > 0 and e=" + e + " < " + ar.length + ".");
		}
	}
	
	/**
	 * checks if the start and end indices are within array bounds
	 * @param ar
	 * @param s
	 * @param e
	 */
	public static void checkForIndicesInBounds(int[] ar, int s, int e) {
		if (s < 0 || e >= ar.length) {
			throw new IndexOutOfBoundsException("Index out of Bounds, s=" + s + " > 0 and e=" + e + " < " + ar.length + ".");
		}
	}
	
	/**
	 * checks if the start and end indices are within array bounds
	 * @param ar
	 * @param s
	 * @param e
	 */
	public static void checkForIndicesInBounds(String[] ar, int s, int e) {
		if (s < 0 || e >= ar.length) {
			throw new IndexOutOfBoundsException("Index out of Bounds, s=" + s + " > 0 and e=" + e + " < " + ar.length + ".");
		}
	}
	
	/**
	 * check if argument not NULL
	 * @param ar
	 */
	public static void checkForNull(double[] ar) {
		if (ar == null) {
			throw new IllegalArgumentException("ar must not be null");
		}
	}
	
	/**
	 * check if arguments are not NULL
	 * @param ar
	 */
	public static void checkForNull2(double[] ... ar) {
		for (int a = 0; a < ar.length; a++) {
			if (ar[a] == null) {
				throw new IllegalArgumentException("arrays must not be NULL");
			}
		}
	}
	
	/**
	 * checks for equal dimensions of both arguments and throws {@code IllegalArgumentException} if not true
	 * @param ar1
	 * @param ar2
	 */
	public static void checkForEqualDimensions(double[] ar1, double[] ar2) {
		if (ar1.length == ar2.length) {
			// do nothing
		} else {
			throw new IllegalArgumentException("ar1[" + ar1.length + "] and ar2[" + ar2.length + "] do not have the same length");
		}
	}
		
	/**
	 * checks for equal dimensions of both arguments and throws {@code IllegalArgumentException} if not true
	 * @param ar1
	 * @param ar2
	 */
	public static void checkForEqualDimensions(float[] ar1, float[] ar2) {
		if (ar1.length == ar2.length) {
			// do nothing
		} else {
			throw new IllegalArgumentException("ar1[" + ar1.length + "] and ar2[" + ar2.length + "] do not have the same length");
		}
	}
	
	/**
	 * checks for equal dimensions of both arguments and throws {@code IllegalArgumentException} if not true
	 * @param ar1
	 * @param ar2
	 */
	public static void checkForEqualDimensions(int[] ar1, int[] ar2) {
		if (ar1.length == ar2.length) {
			// do nothing
		} else {
			throw new IllegalArgumentException("ar1[" + ar1.length + "] and ar2[" + ar2.length + "] do not have the same length");
		}
	}

	/**
	 * checks if the {@code ar} is strictly monotone increasing or decreasing
	 * @param ar
	 * @param increasing
	 * @return
	 */
	public static boolean isStrictMonotone(double[] ar, boolean increasing) {
		if (ar.length < 2) {
			throw new IllegalArgumentException("ar must at least have to elements");
		}
		if (increasing) {
			for (int i = 0; i < ar.length - 1; i++) {
				if (ar[i + 1] <= ar[i]) {
					return false;
				}
			}
		} else {
			for (int i = 0; i < ar.length - 1; i++) {
				if (ar[i + 1] >= ar[i]) {
					return false;
				}
			}
		}
		return true;
	}
	
	/**
	 * checks if the {@code ar} is monotone increasing or decreasing
	 * @param ar
	 * @param increasing
	 * @return
	 */
	public static boolean isMonotone(double[] ar, boolean increasing) {
		if (ar.length < 2) {
			throw new IllegalArgumentException("ar must at least have to elements");
		}
		if (increasing) {
			for (int i = 0; i < ar.length - 1; i++) {
				if (ar[i + 1] < ar[i]) {
					return false;
				}
			}
		} else {
			for (int i = 0; i < ar.length - 1; i++) {
				if (ar[i + 1] > ar[i]) {
					return false;
				}
			}
		}
		return true;
	}
}
