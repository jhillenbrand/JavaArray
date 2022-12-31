package net.sytes.botg.array.math;

import java.util.Arrays;

import net.sytes.botg.array.ArUtils;

public class Vec2Scalar {

	// Suppress default constructor for noninstantiability
	private Vec2Scalar() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * retrieves the maximum in {@code x}
	 * @param x
	 * @return
	 */
	public static double max(double[] x) {
		ArUtils.checkForNull(x);
		ArUtils.checkForEmpty(x);
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
	 * retrieves the maximum valueand the index it was found in {@code x}
	 * <br>both results are stored in double[], where [0] -&gt; max value and [1] -&gt; index
	 * @param x
	 * @return
	 */
	public static double[] maxI(double[] x) {
		ArUtils.checkForNull(x);
		ArUtils.checkForEmpty(x);
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
	 * retrieves the minimum in {@code x}
	 * @param x
	 * @return
	 */
	public static double min(double[] x) {
		ArUtils.checkForNull(x);
		ArUtils.checkForEmpty(x);
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
		ArUtils.checkForEqualDimensions(x1, x2);
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
	public static double scalarProd(double[] x1, double[] x2) {
		ArUtils.checkForEqualDimensions(x1, x2);
		double scalar = 0.0;
		for (int i = 0; i < x1.length; i++) {
			scalar = scalar + x1[i] * x2[i];
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
		ArUtils.checkForNull(p1);
		ArUtils.checkForNull(p2);
		ArUtils.checkForEqualDimensions(p1, p2);
		int n = p1.length;
		double sum = 0;
		for (int i = 0; i < n; i++) {
			sum = sum + Math.pow(p1[i] - p2[i], 2);
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
	 * computes the norm (length) of the vector defined by [x, y]
	 * @param x
	 * @param y
	 * @return
	 */
	public static double norm(double x, double y) {
		return norm(x, y, 0.0);
	}
	
	/**
	 * computes the norm (length) of the vector defined by [x, y, z]
	 * @param x
	 * @param y
	 * @param z
	 * @return
	 */
	public static double norm(double x, double y, double z) {
		return Math.sqrt(x * x + y * y + z * z);
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
		ArUtils.checkForEqualDimensions(x1, x2);
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
		ArUtils.checkForEqualDimensions(x1, x2);
		double sum = 0.0;
		for (int i = 0; i < x1.length; i++) {			
			sum = sum + Math.pow(x1[i] - x2[i], 2);
		}
		return (double) sum / x1.length;
	}
}
