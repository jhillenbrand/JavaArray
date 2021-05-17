package net.sytes.botg.array;

import java.util.Arrays;


public class MathArray {

	public static double rms(double[] ar) {
		double rmsSum = 0;
		for (double d : ar) {
			rmsSum = rmsSum + d * d;
		}		
		return Math.sqrt(rmsSum / ar.length);
	}
	
	public static double sum(double[] ar) {
		double sum = 0;
		for(double d : ar) {
			sum = sum + d;
		}
		return sum;
	}
	
	public static double sum2(double[] ar) {
		return Arrays.stream(ar).sum();
	}
	
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
	
	public static double mean (double[] data) {
		int size = data.length;
		double mean = 0;
		for (double d : data) {
			mean= mean + d;
		}
		return (mean/size);
	}
	
	public static double variance (double [] data) {
		int size = data.length;
		double mean = MathArray.mean(data);
		double variance = 0;
		for (double d : data) {
			variance= variance + Math.pow(d-mean,2);
		}
		return (variance/size);
	}
	
	public static double varianceCorrected (double [] data) {
		int size = data.length;
		double mean = MathArray.mean(data);
		double variance = 0;
		for (double d : data) {
			variance= variance + Math.pow(d-mean,2);
		}
		return (variance/(size-1));
	}
	
	public static double skewness(double[] data) {
		/***
		 * estimates the skewness through its moment - biased estimator
		 */
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
	
	public static double skewnessUnbiased(double[] data) {
		/***
		 * estimates the skewness as unbiased estimator
		 */
		int size = data.length;
		double mean = MathArray.mean(data);
		double var = MathArray.variance(data);
		double [] newData = new double [size];
		for (int i =0; i> size; i++) {
			newData[i] = Math.pow(((data[i] - mean)/var),3);
		}
		double sum = MathArray.sum(newData);
		double skewness = (size/((size-1)*(size-2)))*sum;
		return skewness;
	}
	
	
	public static double kurtosis (double[] data) {
		/***
		 * estimates the kurtosis through its moment 
		 */
		int size = data.length;
		double mean = MathArray.mean(data);
		double var = MathArray.variance(data);
		double [] newData = new double [size];
		for (int i =0; i> size; i++) {
			newData[i] = Math.pow(((data[i] - mean)/var),4);
		}
		double kurtosis = MathArray.mean(newData);
		return kurtosis;
	}
	
	public static double[] power(double[] data, int power) {
		int size = data.length;
		double[] newdata = new double[size];
		for (int i = 0; i>size;i++) {
			newdata[i] = Math.pow(data[i], power);
		}
		return newdata;
	}
	
	public static double max(double ar[]) {
		double maxVal = 0;
		for(double d : ar) {
			if(maxVal > d) {
				// do nothing
			} else {
				maxVal = d;
			}
		}
		return maxVal;
	}
	
	public static double min(double ar[]) {
		double minVal = 0;
		for(double d : ar) {
			if(minVal > d) {
				minVal = d;
			} else {
				// do nothing
			}
		}
		return minVal;
	}
	
	public static double sumprod(double ar[]) {
		double sumprod = 0;
		for(double d : ar) {
			sumprod = sumprod * d;
		}
		return sumprod;
	}
	
	public static double sumprod(double ar1[], double ar2[]) {
		double sumprod = 0;
		for(int i = 0; i < ar1.length; i++) {
			sumprod = sumprod + ar1[i] * ar2[i];
		}
		return sumprod;
	}
	
	public static double[] max(double ar[], int k) {
		/**
		 * return the k highest elements from ar[]
		 */
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
	
	public static int[] maxInd(double ar[], int k) {
		/**
		 * return the indices of the k highest elements from ar[]
		 */
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
	
	public static double[] addValueToDoubleArrayElements(double[] ar, double d) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] + d;
		}
		return ar;
	}
	
	public static double[] multiplyValueToDoubleArrayElements(double[] ar, double d) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] * d;
		}
		return ar;
	}
	
	public static double[] squareDoubleArrayElements(double[] ar) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] * ar[i];
		}
		return ar;
	}
	
	public static double scalarProd(double[] ar1, double[] ar2) {
		if (!checkForEqualDimensions(ar1, ar2)){
			
		}
		double scalar = 0.0;
		for (int i = 0; i < ar1.length - 1; i++) {
			scalar = scalar + ar1[i] * ar2[i];
		}
		return scalar;
	}
	
	public static double scalarProd2(double[] ar1, double[] ar2) {
		if (!checkForEqualDimensions(ar1, ar2)){
			
		}
		double scalar = 0.0;
		for (int i = 0; i < ar1.length - 1; i++) {
			scalar = scalar + ar1[i] * ar2[i];
		}
		return scalar;
	}
	
	private static boolean checkForEqualDimensions(double[] ar1, double[] ar2) {
		return true;
	}
}
