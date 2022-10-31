package net.sytes.botg.array.math;

import net.sytes.botg.array.ArUtils;

public class Mat2Scalar {
	
	// Suppress default constructor for noninstantiability
	private Mat2Scalar() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * return max value of all matrix elements
	 * @param X
	 * @return
	 */
	public static double max(double[][] X) {
		double xMax = Double.MIN_VALUE;
		int n = X.length;
		int m = X[0].length;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (X[i][j] > xMax) {
					xMax = X[i][j];
				}
			}
		}
		return xMax;
	}
	
	/**
	 * return min value of all matrix elements
	 * @param X
	 * @return
	 */
	public static double min(double[][] X) {
		ArUtils.checkForNull(X);
		double xMin = Double.MAX_VALUE;
		int n = X.length;
		int m = X[0].length;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (X[i][j] < xMin) {
					xMin = X[i][j];
				}
			}
		}
		return xMin;
	}
	
	/**
	 * return min and max as array elements of all matrix elements,
	 * <br>Example:
	 * <br>
	 * <br>double[] mm = minMax(x);
	 * <br>min = mm[0];
	 * <br>max = mm[1]; 
	 * @param X
	 * @return
	 */
	public static double[] minMax(double[][] X) {
		ArUtils.checkForNull(X);
		double xMin = Double.MAX_VALUE;
		double xMax = Double.MIN_VALUE;
		int n = X.length;
		int m = X[0].length;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (X[i][j] < xMin) {
					xMin = X[i][j];
				}
				if (X[i][j] > xMax) {
					xMax = X[i][j];
				}
			}
		}
		return new double[]{xMin, xMax};
	}
	
	public static double mean(double[][] X) {
		ArUtils.checkForNull(X);
		double sum = 0.0;
		int n = X.length;
		int m = X[0].length;
		int s = n * m;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				sum = sum + X[i][j];
			}
		}
		return sum / s;
	}
	
	public static double variance(double[][] X) {
		ArUtils.checkForNull(X);
		int n = X.length;
		int m = X[0].length;
		double mean = mean(X);
		double variance = 0.0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				variance = variance + Math.pow(X[i][j] - mean, 2);
			}
		}
		return variance / (n * m);
	}
	
}
