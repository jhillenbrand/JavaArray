package net.sytes.botg.array.math;

import net.sytes.botg.array.ArUtils;

public class Vec2Mat {

	// Suppress default constructor for noninstantiability
	private Vec2Mat() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}	
	
	/**
	 * reshapes the vector into matrix defined by {@code rows}
	 * <br>the elements are sorted column by column
	 * @param x
	 * @param rows
	 * @return
	 */
	public static double[][] vec2Mat(double[] x, int rows){
		return vec2Mat(x, rows, true);
	}
	
	/**
	 * reshapes the vector into matrix defined by {@code rows}
	 * <br>the elements are sorted column by column if {@code columnByColumn} is set to true
	 * @param x
	 * @param rows
	 * @param columnByColumn
	 * @return
	 */
	public static double[][] vec2Mat(double[] x, int rows, boolean columnByColumn){
		ArUtils.checkForNull(x);
		int n = x.length;
		if (n % rows != 0) {
			throw new IllegalArgumentException("the length of vector x (" + n + ") must be divisble by the number of rows (" + rows + ")");
		}
		int columns = n / rows;
		return vec2Mat(x, rows, columns, columnByColumn);
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
	public static double[][] vec2Mat(double[] x, int rows, int columns, boolean columnByColumn){
		ArUtils.checkForNull(x);
		ArUtils.checkForEmpty(x);
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
					double[] nanAr = ArUtils.nan(rem);
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
	
}
