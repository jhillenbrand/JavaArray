package net.sytes.botg.array.math;

import net.sytes.botg.array.ArUtils;

public class Mat2Scalar {
	
	// Suppress default constructor for noninstantiability
	private Mat2Scalar() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * return max value of all matrix elements
	 * @param x
	 * @return
	 */
	public static double max(double[][] x) {
		double xMax = Double.MIN_VALUE;
		int n = x.length;
		int m = x[0].length;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (x[i][j] > xMax) {
					xMax = x[i][j];
				}
			}
		}
		return xMax;
	}
	
	/**
	 * return min value of all matrix elements
	 * @param x
	 * @return
	 */
	public static double min(double[][] x) {
		ArUtils.checkForNull(x);
		double xMin = Double.MAX_VALUE;
		int n = x.length;
		int m = x[0].length;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (x[i][j] < xMin) {
					xMin = x[i][j];
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
	 * @param x
	 * @return
	 */
	public static double[] minMax(double[][] x) {
		ArUtils.checkForNull(x);
		double xMin = Double.MAX_VALUE;
		double xMax = Double.MIN_VALUE;
		int n = x.length;
		int m = x[0].length;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (x[i][j] < xMin) {
					xMin = x[i][j];
				}
				if (x[i][j] > xMax) {
					xMax = x[i][j];
				}
			}
		}
		return new double[]{xMin, xMax};
	}
	
	public static double mean(double[][] x) {
		ArUtils.checkForNull(x);
		double sum = 0.0;
		int n = x.length;
		int m = x[0].length;
		int s = n * m;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				sum = sum + x[i][j];
			}
		}
		return sum / s;
	}
	
	public static double variance(double[][] x) {
		ArUtils.checkForNull(x);
		int n = x.length;
		int m = x[0].length;
		double mean = mean(x);
		double variance = 0.0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				variance = variance + Math.pow(x[i][j] - mean, 2);
			}
		}
		return variance / (n * m);
	}
	
}
