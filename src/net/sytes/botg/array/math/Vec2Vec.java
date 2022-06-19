package net.sytes.botg.array.math;

import java.util.Arrays;

import net.sytes.botg.array.ArrayUtility;
import net.sytes.botg.array.SearchArray;
import net.sytes.botg.array.SortArray;

public class Vec2Vec {

	// Suppress default constructor for noninstantiability
	private Vec2Vec() {
		throw new AssertionError();
	}
		

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
		return offset(data, -Vec2Scalar.mean(data));
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
		ArrayUtility.checkForFirstSmallerSecond(lowerLimit, upperLimit);
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
		ArrayUtility.checkForFirstSmallerSecond(lowerLimit, upperLimit);
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
		ArrayUtility.checkForEqualDimensions(t, ar);
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
		ArrayUtility.checkForEqualDimensions(ar1, ar2);
		int L = ar1.length;
		double[] ar3 = new double[L];
		for (int i = 0; i < L; i++) {
			ar3[i] = ar1[i] + ar2[i];
		}
		return ar3;
	}
	
	/**
	 * element wise subtraction of array values
	 * @param ar1
	 * @param ar2
	 * @return
	 */
	public static double[] minus(double[] ar1, double[] ar2){
		ArrayUtility.checkForEqualDimensions(ar1, ar2);
		int L = ar1.length;
		double[] ar3 = new double[L];
		for (int i = 0; i < L; i++) {
			ar3[i] = ar1[i] - ar2[i];
		}
		return ar3;
	}

	/**
	 * add the element x to the end of the array ar and return the new array
	 * @param <T>
	 * @param x
	 * @param ar
	 * @return
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
}
