package net.sytes.botg.array.math;

import java.util.Arrays;

import net.sytes.botg.array.SearchArray;
import net.sytes.botg.array.SortArray;

public class MathArray {

	/**
	 * ---------------------------------------------------------------------------
	 * SECTION - WINDOWING
	 */
	
	/**
	 * separates the given double array into windows of size windowSize
	 * <br>if dropUneven is true, then the remaining last ar.length % windowSize elements are dropped
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
	 * ENDSECTION - WINDOWING
	 * ----------------------------------------------------------------------------
	 */
	
	
	/**
	 * ---------------------------------------------------------------------------
	 * SECTION - SIGNAL CORRECTIONS
	 */
	
	/**
	 * methods removes the offset and mirrors negative values above ZERO
	 * @param data
	 * @return
	 */
	public static double[] rectify(double[] data) {
		return abs(removeOffset(data));
	}
	
	/**
	 * removes the offset of {@code data}
	 * @param data
	 * @return
	 */
	public static double[] removeOffset(double[] data) {
		return offset(data, -mean(data));
	}
	
	/**
	 * removes {@code NaN} entries in {@code data}
	* @param data
	 * @return double[]
	 */
	public static double[] removeNaN(double[] data) {
		return SearchArray.elementsAt(SearchArray.isNaN(data), data);
	}
	
	/**
	 * returns a double[] array containing all elements in {@code data} within the bounds [{@code lowerLimit}, {@code upperLimit}]
	 * @param data
	 * @param lowerLimit
	 * @param upperLimit
	 * @return double[]
	 */
	public static double[] removeDataOutOfRange(double[] data, double lowerLimit, double upperLimit) {
		checkForFirstSmallerSecond(lowerLimit, upperLimit);
		boolean[] inRange = SearchArray.isInRange(data, lowerLimit, upperLimit);
		int n = sum(inRange);
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
		checkForFirstSmallerSecond(lowerLimit, upperLimit);
		boolean[] inRange = SearchArray.isInRange(data, lowerLimit, upperLimit);
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
	 * ENDSECTION - SIGNAL CORRECTIONS
	 * ----------------------------------------------------------------------------
	 */
	
	/**
	 * If a number is divisible by 2 then it has its least significant bit (LSB) set to 0,
	 * 	if divisible by 4 then two LSB’s set to 0, if by 8 then three LSB’s set to 0 and so on.
	 * 	Keeping this in mind, a number n is divisible by 2m if (n & ((1 << m) – 1)) is equal to 0 else not.
	 * @param n
	 * @param m
	 * @return
	 */
	public static boolean isDivByPow(int n, int m) {
		if ((n & ((1 << m) - 1)) == 0) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * computes the closest exponent for base 2 for {@code n}
	 * @param n
	 * @return
	 */
	public static int closestExponentForBase2(int n) {
		int m = (int) log2(n);		
		return m;
	}
	
	/**
	 * computes the log of base 2 of @code d
	 * @param d
	 * @return
	 */
	public static double log2(double d) {
		return Math.log(d) / Math.log(2);
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
	public static double sum(int[] ar) {
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
	public static double variance (double [] data) {
		int size = data.length;
		double mean = MathArray.mean(data);
		double variance = 0;
		for (double d : data) {
			variance= variance + Math.pow(d-mean,2);
		}
		return (variance/size);
	}
	
	/**
	 * estimates the corrected variance
	 * @param data
	 * @return
	 */
	public static double varianceCorrected (double [] data) {
		int size = data.length;
		double mean = MathArray.mean(data);
		double variance = 0;
		for (double d : data) {
			variance= variance + Math.pow(d-mean,2);
		}
		return (variance/(size-1));
	}
	
	/**
	 * estimates the skewness through its moment - biased estimator
	 * @param data
	 * @return
	 */
	public static double skewness(double[] data) {
		int size = data.length;
		double mean = MathArray.mean(data);
		double var = MathArray.variance(data);
		double [] newData = new double [size];
		for (int i =0; i> size; i++) {
			newData[i] = Math.pow(((data[i] - mean)/var),3);
		}
		double skewness = MathArray.mean(newData);
		return skewness;
	}
	
	/**
	 * estimates the skewness as unbiased estimator
	 * @param ar
	 * @return
	 */
	public static double skewnessUnbiased(double[] ar) {
		int size = ar.length;
		double mean = MathArray.mean(ar);
		double var = MathArray.variance(ar);
		double [] newAr = new double [size];
		for (int i =0; i> size; i++) {
			newAr[i] = Math.pow(((ar[i] - mean)/var),3);
		}
		double sum = MathArray.sum(newAr);
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
		double mean = MathArray.mean(ar);
		double var = MathArray.variance(ar);
		double [] newAr = new double [size];
		for (int i =0; i> size; i++) {
			newAr[i] = Math.pow(((ar[i] - mean)/var),4);
		}
		double kurtosis = MathArray.mean(newAr);
		return kurtosis;
	}
	
	/**
	 * takes every element in {@code ar} to the power of {@code p}
	 * @param data
	 * @param power
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
	 * retrieves the maximum in {@code ar}
	 * @param ar
	 * @return
	 */
	public static double max(double[] ar) {
		checkForNull(ar);
		checkForEmpty(ar);
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
	 * retrieves the minimum in {@code ar}
	 * @param ar
	 * @return
	 */
	public static double min(double[] ar) {
		checkForNull(ar);
		checkForEmpty(ar);
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
	 *<br>span({1.0, 2.5, 3.4}) --> 3.4 - 1.0 = 2.4 
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
	 * computes the sum product of {@code ar1} with {@code ar2}
	 * @param ar1
	 * @param ar2
	 * @return
	 */
	public static double sumprod(double[] ar1, double[] ar2) {
		checkForEqualDimensions(ar1, ar2);
		double sumprod = 0;
		for(int i = 0; i < ar1.length; i++) {
			sumprod = sumprod + ar1[i] * ar2[i];
		}
		return sumprod;
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
	 * return the indices of the k highest elements from ar[]
	 * @param ar
	 * @param k
	 * @return
	 */
	public static int[] maxInd(double[] ar, int k) {
		int[] sortInds = SortArray.quicksort2(ar);
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
	
	/**
	 * offsets every element in {@code ar} by {@code o}  
	 * @param ar
	 * @param o
	 * @return double[]
	 */
	public static double[] offset(double[] ar, double o) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] + o;
		}
		return ar;
	}
	
	/**
	 * offsets every element in {@code ar} by {@code o} 
	 * @param ar
	 * @param o
	 * @return
	 */
	public static long[] offset(long[] ar, long o) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] + o;
		}
		return ar;
	}
	
	/**
	 * offsets every element in {@code ar} by {@code o} 
	 * @param ar
	 * @param o
	 * @return
	 */
	public static int[] offset(int[] ar, int o) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] + o;
		}
		return ar;
	}
	
	/**
	 * multiplies {@code d} to every element in {@code ar}
	 * @param ar
	 * @param d
	 * @return
	 */
	public static double[] multiply(double[] ar, double d) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] * d;
		}
		return ar;
	}
	
	/**
	 * squares every element in {@code ar}
	 * @param ar
	 * @return
	 */
	public static double[] square(double[] ar) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] * ar[i];
		}
		return ar;
	}
	
	/**
	 * computes the scalar product between {@code ar1} and {@code ar2}
	 * @param ar1
	 * @param ar2
	 * @return
	 */
	public static double scalarProd(double[] ar1, double[] ar2) {
		checkForEqualDimensions(ar1, ar2);
		double scalar = 0.0;
		for (int i = 0; i < ar1.length; i++) {
			scalar = scalar + ar1[i] * ar2[i];
		}
		return scalar;
	}
	
	/**
	 * checks for equal dimensions of both arguments and throws {@code IllegalArgumentException} if not true
	 * @param ar1
	 * @param ar2
	 */
	private static void checkForEqualDimensions(double[] ar1, double[] ar2) {
		if (ar1.length == ar2.length) {
			// do nothing
		} else {
			throw new IllegalArgumentException("ar1[" + ar1.length + "] and ar2[" + ar2.length + "] do not have the same length");
		}
	}
	
	/**
	 * checks if the first argument is smaller than the second and throws {@code IllegalArgumentException} if true
	 * @param d1
	 * @param d2
	 */
	private static void checkForFirstSmallerSecond(double d1, double d2) {
		if (d1 >= d2) {
			throw new IllegalArgumentException("d1[" + d1 + "] must be smaller than d2[" + d2 + "]");
		}
	}
	
	private static void checkForNull(double[] ar) {
		if (ar == null) {
			throw new IllegalArgumentException("ar must not be null");
		}
	}
	
	private static void checkForEmpty(double[] ar) {
		if (ar.length == 0) {
			throw new IllegalArgumentException("ar must not be empty");
		}
	}
}
