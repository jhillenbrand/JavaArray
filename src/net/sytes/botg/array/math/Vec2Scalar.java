package net.sytes.botg.array.math;

import java.util.Arrays;

import net.sytes.botg.array.ArUtils;

public class Vec2Scalar {

	// Suppress default constructor for noninstantiability
	private Vec2Scalar() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * retrieves the maximum in {@code ar}
	 * @param ar
	 * @return
	 */
	public static double max(double[] ar) {
		ArUtils.checkForNull(ar);
		ArUtils.checkForEmpty(ar);
		double maxVal = ar[0];
		for(double d : ar) {
			if(maxVal > d) {
				// do nothing
			} else {
				maxVal = d;
			}
		}
		return maxVal;
	}
	
	/**
	 * retrieves the maximum valueand the index it was found in {@code ar}
	 * <br>both results are stored in double[], where [0] -&gt; max value and [1] -&gt; index
	 * @param ar
	 * @return
	 */
	public static double[] maxI(double[] ar) {
		ArUtils.checkForNull(ar);
		ArUtils.checkForEmpty(ar);
		double maxVal = ar[0];
		int maxInd = -1;
		for(int i = 0; i < ar.length; i++) {
			if(maxVal > ar[i]) {
				// do nothing
			} else {
				maxVal = ar[i];
				maxInd = i;
			}
		}
		return new double[] {maxVal, maxInd};
	}
	
	/**
	 * retrieves the minimum in {@code ar}
	 * @param ar
	 * @return
	 */
	public static double min(double[] ar) {
		ArUtils.checkForNull(ar);
		ArUtils.checkForEmpty(ar);
		double minVal = ar[0];
		for(double d : ar) {
			if(minVal > d) {
				minVal = d;
			} else {
				// do nothing
			}
		}
		return minVal;
	}
	
	/**
	 * returns the span value of {@code ar}
	 *<br>Example:
	 *<br>span({1.0, 2.5, 3.4}) --&gt; 3.4 - 1.0 = 2.4 
	 * @param ar
	 * @return
	 */
	public static double span(double[] ar) {
		double maxVal = ar[0];
		double minVal = ar[0];
		for(double d : ar) {
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
	 * rms value of {@code ar}
	 * @param ar
	 * @return
	 */
	public static double rms(double[] ar) {
		double rmsSum = 0;
		for (double d : ar) {
			rmsSum = rmsSum + d * d;
		}		
		return Math.sqrt(rmsSum / ar.length);
	}
	
	/**
	 * rms mean value of {@code ar}
	 * @param ar
	 * @return
	 */
	public static double rmsMean(double[] ar) {
		double rmsSum = 0;
		double sum = 0;
		for (double d : ar) {
			rmsSum = rmsSum + d * d;
			sum = sum + d;
		}		
		return Math.sqrt(rmsSum / ar.length) + sum / ar.length;
	}
	
	/**
	 * sum of @code ar	 * 
	 * <br>is way faster than sum2
	 * @param ar
	 * @return
	 */
	public static double sum(double[] ar) {
		double sum = 0;
		for(double d : ar) {
			sum = sum + d;
		}
		return sum;
	}
	
	/**
	 * returns the sum of Object primitive Double[]
	 * @param ar
	 * @return
	 */
	public static double sum(Double[] ar) {
		double sum = 0;
		for(double d : ar) {
			sum = sum + d;
		}
		return sum;
	}
	
	/**
	 * sum of @code ar
	 * @param ar
	 * @return
	 */
	public static double sum (int[] ar) {
		int sum = 0;
		for(int d : ar) {
			sum = sum + d;
		}
		return sum;
	}
	
	/**
	 * sum of @code ar
	 * @param ar
	 * @return
	 */
	public static long sum(long[] ar) {
		long sum = 0;
		for(long lo : ar) {
			sum = sum + lo;
		}
		return sum;
	}
	
	/**
	 * sum of {@code ar} based on Streams
	 * @param ar
	 * @return
	 */
	public static double sum2(double[] ar) {
		return Arrays.stream(ar).sum();
	}
	
	/**
	 * counts the number of {@code} true elements in {@code bools}
	 * @param bools
	 * @return int
	 */
	public static int sum(boolean[] bools) {
		int sum = 0;
		for (boolean b : bools) {
		    sum += b ? 1 : 0;
		}
		return sum;
	}
	
	/**
	 * computes the median
	 * @param ar
	 * @return
	 */
	public static double median(double ar[]) {
		Arrays.sort(ar);
		double med = 0;
		if (ar.length % 2 == 0) {
		    med = (ar[ar.length / 2] + ar[ar.length / 2 - 1]) / 2;
		} else {
		    med = ar[ar.length / 2];
		}
		return med;
	}
	
	/**
	 * computes the mean value of the given {@code data}
	 * @param data
	 * @return double
	 */
	public static double mean(double[] data) {
		return sum(data) / data.length;
	}
	
	/**
	 * computes the mean value of the given {@code data} as Object Primitives Double[]
	 * @param data
	 * @return
	 */
	public static double mean(Double[] data) {
		return sum(data) / data.length;
	}
	
	/**
	 * computes the mean of long Array
	 * @param data
	 * @return
	 */
	public static double mean(long[] data) {
		return (double) sum(data) / data.length;
	}
	
	/**
	 * estimates the variance
	 * @param data
	 * @return
	 */
	public static double variance (double[] data) {
		int size = data.length;
		double mean = mean(data);
		double variance = 0;
		for (double d : data) {
			variance = variance + Math.pow(d - mean, 2);
		}
		return variance / size;
	}
	
	/**
	 * estimates the corrected variance
	 * @param data
	 * @return
	 */
	public static double varianceCorrected (double [] data) {
		int size = data.length;
		double mean = mean(data);
		double variance = 0;
		for (double d : data) {
			variance = variance + Math.pow(d - mean, 2);
		}
		return variance / (size - 1);
	}
	
	/**
	 * estimates the skewness through its moment - biased estimator
	 * @param data
	 * @return
	 */
	public static double skewness(double[] data) {
		int size = data.length;
		double mean = mean(data);
		double var = variance(data);
		double [] newData = new double [size];
		for (int i =0; i> size; i++) {
			newData[i] = Math.pow(((data[i] - mean) / var), 3);
		}
		double skewness = mean(newData);
		return skewness;
	}
	
	/**
	 * estimates the skewness as unbiased estimator
	 * @param ar
	 * @return
	 */
	public static double skewnessUnbiased(double[] ar) {
		int size = ar.length;
		double mean = mean(ar);
		double var = variance(ar);
		double [] newAr = new double [size];
		for (int i =0; i> size; i++) {
			newAr[i] = Math.pow(((ar[i] - mean)/var),3);
		}
		double sum = sum(newAr);
		double skewness = (size/((size-1)*(size-2)))*sum;
		return skewness;
	}
	
	/**
	 * estimates the kurtosis through its moment 
	 * @param ar
	 * @return
	 */
	public static double kurtosis (double[] ar) {
		int size = ar.length;
		double mean = mean(ar);
		double var = variance(ar);
		double [] newAr = new double [size];
		for (int i =0; i> size; i++) {
			newAr[i] = Math.pow(((ar[i] - mean)/var),4);
		}
		double kurtosis = mean(newAr);
		return kurtosis;
	}
		
	/**
	 * computes the sum product of {@code ar1} with {@code ar2}
	 * @param ar1
	 * @param ar2
	 * @return
	 */
	public static double sumprod(double[] ar1, double[] ar2) {
		ArUtils.checkForEqualDimensions(ar1, ar2);
		double sumprod = 0;
		for(int i = 0; i < ar1.length; i++) {
			sumprod = sumprod + ar1[i] * ar2[i];
		}
		return sumprod;
	}
	
	/**
	 * computes the scalar product between {@code ar1} and {@code ar2}
	 * @param ar1
	 * @param ar2
	 * @return
	 */
	public static double scalarProd(double[] ar1, double[] ar2) {
		ArUtils.checkForEqualDimensions(ar1, ar2);
		double scalar = 0.0;
		for (int i = 0; i < ar1.length; i++) {
			scalar = scalar + ar1[i] * ar2[i];
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
	 * computes the mean squared error between two double arrays
	 * @param ar1
	 * @param ar2
	 * @return
	 */
	public static double mse(double[] ar1, double[] ar2) {
		ArUtils.checkForEqualDimensions(ar1, ar2);
		double sum = 0.0;
		for (int i = 0; i < ar1.length; i++) {			
			sum = sum + Math.pow(ar1[i] - ar2[i], 2);
		}
		return sum / ar1.length;
	}
	
	/**
	 * computes the mean squared error between two int arrays
	 * @param ar1
	 * @param ar2
	 * @return 
	 */
	public static double mse(int[] ar1, int[] ar2) {
		ArUtils.checkForEqualDimensions(ar1, ar2);
		double sum = 0.0;
		for (int i = 0; i < ar1.length; i++) {			
			sum = sum + Math.pow(ar1[i] - ar2[i], 2);
		}
		return (double) sum / ar1.length;
	}
}
