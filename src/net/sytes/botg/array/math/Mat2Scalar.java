package net.sytes.botg.array.math;

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
		double xMin = Double.MAX_VALUE;
		int n = x.length;
		int m = x[0].length;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (x[i][j] > xMin) {
					xMin = x[i][j];
				}
			}
		}
		return xMin;
	}
	
	public static double mean(double[][] x) {
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
		
	}
	
}
