package net.sytes.botg.array.math;

import net.sytes.botg.array.ArUtils;

public class Mat2Vec {

	// Suppress default constructor for noninstantiability
	private Mat2Vec() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * returns the matrix {@code X} as vector, where all matrix elements are listed column by column
	 * @param X
	 * @return
	 */
	public static double[] mat2Vec(double[][] X) {
		return mat2Vec(X, true);
	}
	
	/**
	 * 
	 * @param X
	 * @return
	 */
	public static double[] mat2Vec(double[][] X, boolean columnByColumn) {
		ArUtils.checkForNull(X);
		ArUtils.checkForEmpty(X);
		int n = X.length;
		int m = X[0].length;
		double[] x = new double[n * m];
		if (columnByColumn) {
			int c = 0;
			for (int j = 0; j < m; j++) {
				for (int i = 0; i < n; i++) {
					x[c] = X[i][j];
					++c;
				}
			}
		} else {
			int c = 0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					x[c] = X[i][j];
					++c;
				}
			}
		}
		return x;
	}
	
}
