package net.sytes.botg.array;

import java.util.Arrays;
import java.util.Random;

/**
 * This class contains helper methods for array creation, function checks for size, indices, etc.
 * <br>and 2D console print methods for MATLAB like output
 * @author hillenbrand
 *
 */
public class ArUtils {

	// Suppress default constructor for noninstantiability
	private ArUtils() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	private static final long BETTER_OF_AS_STREAM_SIZE = 100_000_000;
		
	/**
	 * returns an array of size n with 0's
	 * @param n
	 * @return
	 */
	public static double[] zeros(int n) {
		return new double[n];
	}
	
	/**
	 * returns an float array of size n with 0's
	 * @param n
	 * @return
	 */
	public static float[] zerosF(int n) {
		return new float[n];
	}
	
	/**
	 * return an array with zeros of size [n, m]
	 * @param n
	 * @param m
	 * @return
	 */
	public static double[][] zeros(int n, int m) {
		return new double[n][m];
	}
	
	/**
	 * return an float array with zeros of size [n, m]
	 * @param n
	 * @param m
	 * @return
	 */
	public static float[][] zerosF(int n, int m) {
		return new float[n][m];
	}
	
	/**
	 * returns an array of size {@code n} with 1's
	 * @param n
	 * @return
	 */
	public static double[] ones(int n) {
		double[] ar = new double[n];
		for (int i = 0; i < n; i++) {
			ar[i] = 1;
		}
		return ar;
	}
	
	/**
	 * returns an float array of size {@code n} with 1's
	 * @param n
	 * @return
	 */
	public static float[] onesF(int n) {
		float[] x = new float[n];
		for (int i = 0; i < n; i++) {
			x[i] = 1;
		}
		return x;
	}
	
	/**
	 * returns an int array of size {@code n} with 1's
	 * @param n
	 * @return
	 */
	public static int[] onesI(int n) {
		int[] x = new int[n];
		for (int i = 0; i < n; i++) {
			x[i] = 1;
		}
		return x;
	}
	
	/**
	 * returns an 2D array with 1's, where {@code n} is the number of rows and {@code m} the number of columns
	 * @param n
	 * @param m
	 * @return
	 */
	public static double[][] ones(int n, int m){
		double[][] ar = new double[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				ar[i][j] = 1.0;
			}
		}
		return ar;
	}
	
	/**
	 * returns an 2D float array with 1's, where {@code n} is the number of rows and {@code m} the number of columns
	 * @param n
	 * @param m
	 * @return
	 */
	public static float[][] onesF(int n, int m){
		float[][] ar = new float[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				ar[i][j] = 1.0f;
			}
		}
		return ar;
	}
	
	/**
	 * returns a squared unit matrix of size {@code n}
	 * <br> Example for n = 3
	 * <br>| 1, 0, 0 |
	 * <br>| 0, 1, 0 | 
	 * <br>| 0, 0, 1 |
	 * @param n
	 * @return
	 */
	public static double[][] unitMatrix(int n){
		double[][] ar = zeros(n, n);
		for (int i = 0; i < n; i++) {
			ar[i][i] = 1;
		}
		return ar;
	}
	
	/**
	 * returns a unit square float matrix of size {@code n}
	 * <br> Example for n = 3
	 * <br>| 1f, 0,  0 |
	 * <br>| 0, 1f,  0 | 
	 * <br>| 0,  0, 1f |
	 * @param n
	 * @return
	 */
	public static float[][] unitMatrixF(int n){
		float[][] ar = zerosF(n, n);
		for (int i = 0; i < n; i++) {
			ar[i][i] = 1;
		}
		return ar;
	}
	
	/**
	 * creates a new random double array with values in [0, 1]
	 * @param n
	 * @return
	 */
	public static double[] rand(int n) {
		Random r = new Random();
		double[] ar = new double[n];
		for (int i = 0; i < n; i++) {
			ar[i] = r.nextDouble();
		}
		return ar;
	}
	
	/**
	 * creates a new random double array with values in [{@code start},{@code end}]
	 * @param n
	 * @param start
	 * @param end
	 * @return
	 */
	public static double[] rand(int n, double start, double end) {
		ArUtils.checkForFirstSmallerSecond(start, end);
		Random r = new Random();
		double[] ar = new double[n];
		for (int i = 0; i < n; i++) {
			ar[i] = (end - start) * r.nextDouble() + start;
		}
		return ar;
	}
	
	/**
	 * returns a 2D matrix with random values in [0, 1] of size [{@code n},{@code m}]
	 * @param n
	 * @param m
	 * @return
	 */
	public static double[][] rand(int n, int m){
		Random r = new Random();
		double[][] ar = new double[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				ar[i][j] = r.nextDouble();
			}
		}
		return ar;
	}
	
	public static int randInt(int min, int max) {
		ArUtils.checkForFirstSmallerSecond(min, max);
		Random r = new Random();
		return r.nextInt((max - min) + 1) + min;
	}
	
	public static int[] randInt(int n, int min, int max) {
		if (min >= max) {
			throw new IllegalArgumentException("max must be greater than min");
		}
		int[] rInts = new int[n];		
		Random r = new Random();
		for (int i = 0; i < n; i++) {
			rInts[i] = r.nextInt((max - min) + 1) + min;
		}
		return rInts;
	}
	
	public static int[] randInt(int size) {
		int[] data;
		if (size < BETTER_OF_AS_STREAM_SIZE) {
			data = new int[size];
			Random r = new Random();
			for (int i = 0; i < size; i++) {
				data[i] = r.nextInt();
			}
		} else {
			data = new Random().ints(size).toArray();
		}
		return data;
	}
	
	/**
	 * returns a matrix X € R<sup>nxm</sup> with elements increasing from 0 to k = {@code n} * {@code m}
	 * @param n
	 * @param m
	 * @return
	 */
	public static double[][] incrementMat(int n, int m) {
		int k = n * m;
		double[][] X = new double[n][m];		
		for (int i = 0; i < k; i++) {
			int r = i / m;
			int c = i % m;
			X[r][c] = i;
		}
		return X;
	}
	
	/**
	 * returns a double array starting from {@code start} to {@code end} with {@code step}
	 * @param start
	 * @param end
	 * @param step
	 * @return
	 */
	public static double[] linspace(double start, double end, double step) {
		int size = (int) ((end - start ) / step);
		double[] ar = new double[size];
		for (int i = 0; i < size; i++) {
			ar[i] = start + step * i;
		}
		return ar;
	}
	
	/**
	 * returns a int array starting from {@code start} to {@code end} with {@code step}
	 * @param start
	 * @param end
	 * @param step
	 * @return
	 */
	public static int[] linspace(int start, int end, int step) {
		int size = (int) ((end - start ) / step);
		int[] ar = new int[size];
		for (int i = 0; i < size; i++) {
			ar[i] = start + step * i;
		}
		return ar;
	}
	
	/**
	 * returns a double array starting from {@code start} to {@code end} with equally spaced steps, so that the {@code size} is met
	 * @param start
	 * @param end
	 * @param size
	 * @return
	 */
	public static double[] linspace(double start, double end, int size) {
		double[] ar = new double[size];
		double step = (end - start) / (size - 1);
		for (int i = 0; i < size; i++) {
			ar[i] = start + step * i;
		}
		return ar;
	}
	
	/**
	 * returns a int array starting from {@code start} to {@code end} with step size 1
	 * @param start
	 * @param end
	 * @return
	 */
	public static int[] linspace(int start, int end) {
		int size = (int) ((end - start ) / 1);
		int[] ar = new int[size];
		for (int i = 0; i < size; i++) {
			ar[i] = start + 1 * i;
		}
		return ar;
	}
	
	/**
	 * returns a double array starting from {@code start} with {@code size} steps of {@code step}
	 * @param start
	 * @param size
	 * @param step
	 * @return
	 */
	public static double[] linspace(double start, int size, double step) {
		double[] ar = new double[size];
		ar[0] = start;
		for (int i = 1; i < size; i++) {
			ar[i] = ar[i - 1] + step;
		}
		return ar;
	}
	
	/**
	 * returns a long array starting from {@code start} with size {@code size} and steps of {@code step}
	 * @param start
	 * @param size
	 * @param step
	 * @return
	 */
	public static long[] linspace(long start, int size, long step) {
		long[] ar = new long[size];
		ar[0] = start;
		for (int i = 1; i < size; i++) {
			ar[i] = ar[i - 1] + step;
		}
		return ar;
	}
	
	/**
	 * returns a double array starting from 0 and incrementing by 1 for {@code size} times
	 * @param size
	 * @return
	 */
	public static double[] linspace(int size) {
		return linspace(0, size, 1.0);
	}
	
	/**
	 * Create vector of {@code n} logarithmically spaced values between {@code start} and {@code end} 
	 * @param start
	 * @param end
	 * @param n
	 * @return
	 */
	public static double[] logspace(double start, double end, int n) {
		double[] ar = new double[n];
		double step = (end - start) / (n - 1);
		for (int i = 0; i < n; i++) {
			ar[i] = Math.pow(start + step * i, 10);
		}
		return ar;
	}
	
	/**
	 * creates a meshgrid of the two input vectors
	 * @param x
	 * @param y
	 * @return
	 */
	public static double[][][] meshgrid(double[] x, double[] y){
		
		int n = x.length;
		int m = y.length;
		
		double[][] X = new double[n][m];
		double[][] Y = new double[n][m];
		
		for (int i = 0; i < n; i++) {			
			for (int j = 0; j < m; j++) {
				X[i][j] = x[i];
				Y[i][j] = y[j];
			}			
		}
		
		double[][][] XY = new double[2][][];
		XY[0] = X;
		XY[1] = Y;
		
		return XY;
		
	}
	
	/**
	 * returns a vector with NaN values of length {@code n}
	 * @param size
	 * @return
	 */
	public static double[] nan(int n) {
		double[] data = new double[n];
		for (int i = 0; i < n; i++) {
			data[i] = Double.NaN;
		}
		return data;
	}
	
	/**
	 * returns 2D matrix with NaN values of length {@code n}
	 * @param size
	 * @return
	 */
	public static double[][] nan(int n, int m) {
		double[][] mat = new double[m][n];
		for (int j = 0; j < m; j++) {
			mat[j] = nan(n);
		}
		return mat;
	}
	
	/**
	 * returns a sub array starting at 0 and ending at index e from ar
	 * @param ar
	 * @param e
	 * @return
	 */
	public static double[] subArray(double[] ar, int e) {
		return subArray(ar, 0, e);
	}
	
	/**
	 * returns a sub array starting at index s and ending at index e from ar
	 * @param ar
	 * @param s
	 * @param e
	 * @return
	 */
	public static double[] subArray(final double[] ar, int s, int e) {
		checkForIndicesInBounds(ar, s, e);
		double[] ar2 = new double[e - s + 1];
		System.arraycopy(ar, s, ar2, 0, ar2.length);
		return ar2;
	}
	
	/**
	 * returns a sub array starting at 0 and ending at index e from ar
	 * @param ar
	 * @param e
	 * @return
	 */
	public static int[] subArray(int[] ar, int e) {
		return subArray(ar, 0, e);
	}
	
	/**
	 * returns a sub array starting at index s and ending at index e from ar
	 * @param ar
	 * @param s
	 * @param e
	 * @return
	 */
	public static int[] subArray(final int[] ar, int s, int e) {
		checkForIndicesInBounds(ar, s, e);
		int[] ar2 = new int[e - s + 1];
		System.arraycopy(ar, s, ar2, 0, ar2.length);
		return ar2;
	}
	
	/**
	 * returns a sub array starting at index s and ending at index e from ar
	 * @param ar
	 * @param s
	 * @param e
	 * @return
	 */
	public static String[] subArray(final String[] ar, int s, int e) {
		checkForIndicesInBounds(ar, s, e);
		String[] ar2 = new String[e - s + 1];
		System.arraycopy(ar, s, ar2, 0, ar2.length);
		return ar2;
	}
	
	public static double[] copy(double[] x) {
		return x.clone();
	}
	
	/**
	 * deep clone an matrix {@code X}
	 * @param X
	 * @return
	 */
	public static double[][] copy(double[][] X){
		checkForNull(X);
		checkForEmpty(X);
		int m = X.length;
		int n = X[0].length;
		double[][] Y = new double[m][n];
		for (int j = 0; j < m; j++) {
			System.arraycopy(X[j], 0, Y[j], 0, X[j].length);
		}
		return Y;
	}
	
	/**
	 * deep clone an matrix {@code X} using loop
	 * @param X
	 * @return
	 */
	public static double[][] copy2(double[][] X){
		checkForNull(X);
		checkForEmpty(X);
		int m = X.length;
		int n = X[0].length;
		double[][] Y = new double[m][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				Y[j][i] = X[j][i];
			}
		}
		return Y;
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
	 * checks for equal dimensions of both arguments and throws {@code IllegalArgumentException} if not true
	 * @param ar1
	 * @param ar2
	 */
	public static void checkForEqualDimensions(double[] ar1, double[] ar2) {
		if (ar1.length == ar2.length) {
			// do nothing
		} else {
			throw new IllegalArgumentException("ar1[" + ar1.length + "] and ar2[" + ar2.length + "] do not have the same length");
		}
	}
	
	/**
	 * checks for equal dimensions of both arguments and throws {@code IllegalArgumentException} if not true
	 * @param ar1
	 * @param ar2
	 */
	public static void checkForEqualDimensions(float[] ar1, float[] ar2) {
		if (ar1.length == ar2.length) {
			// do nothing
		} else {
			throw new IllegalArgumentException("ar1[" + ar1.length + "] and ar2[" + ar2.length + "] do not have the same length");
		}
	}
	
	/**
	 * checks for equal dimensions of both arguments and throws {@code IllegalArgumentException} if not true
	 * @param ar1
	 * @param ar2
	 */
	public static void checkForEqualDimensions(int[] ar1, int[] ar2) {
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
	public static void checkForFirstSmallerSecond(double d1, double d2) {
		if (d1 >= d2) {
			throw new IllegalArgumentException("d1[" + d1 + "] must be smaller than d2[" + d2 + "]");
		}
	}
	
	/**
	 * check if argument not NULL
	 * @param ar
	 */
	public static void checkForNull(double[] ar) {
		if (ar == null) {
			throw new IllegalArgumentException("ar must not be null");
		}
	}
	
	/**
	 * check if arguments are not NULL
	 * @param ar
	 */
	public static void checkForNull2(double[] ... ar) {
		for (int a = 0; a < ar.length; a++) {
			if (ar[a] == null) {
				throw new IllegalArgumentException("arrays must not be NULL");
			}
		}
	}
	
	public static void checkForNull(double[][] ar) {
		if (ar == null) {
			throw new IllegalArgumentException("ar must not be null");
		}
	}
	
	public static void checkForNull(double[][] ... ar) {
		for (int i = 0; i < ar.length; i++) {
			if (ar[i] == null) {
				throw new IllegalArgumentException("ar must not be null");
			}
		}
	}
	
	public static void checkForEmpty(double[] ar) {
		if (ar.length == 0) {
			throw new IllegalArgumentException("ar must not be empty");
		}
	}
	
	public static void checkForEmpty(double[][] X) {
		if (X.length == 0) {
			throw new IllegalArgumentException("matrix X must not be empty");
		} else {
			for (int i = 0; i < X.length; i++) {
				if (X[i].length == 0) {
					throw new IllegalArgumentException("matrix X must not contain empty vectors");
				}
			}
		}
	}
	
	public static void checkForEmpty(double[][] ... ar) {
		for (int i = 0; i < ar.length; i++) {
			if (ar.length == 0) {
				throw new IllegalArgumentException("ar must not be empty");
			}
		}
	}
	
	/**
	 * checks if the start and end indices are within array bounds
	 * @param ar
	 * @param s
	 * @param e
	 */
	public static void checkForIndicesInBounds(double[] ar, int s, int e) {
		if (s < 0 || e >= ar.length) {
			throw new IndexOutOfBoundsException("Index out of Bounds, s=" + s + " > 0 and e=" + e + " < " + ar.length + ".");
		}
	}
	
	/**
	 * checks if the start and end indices are within array bounds
	 * @param ar
	 * @param s
	 * @param e
	 */
	public static void checkForIndicesInBounds(int[] ar, int s, int e) {
		if (s < 0 || e >= ar.length) {
			throw new IndexOutOfBoundsException("Index out of Bounds, s=" + s + " > 0 and e=" + e + " < " + ar.length + ".");
		}
	}
	
	/**
	 * checks if the start and end indices are within array bounds
	 * @param ar
	 * @param s
	 * @param e
	 */
	public static void checkForIndicesInBounds(String[] ar, int s, int e) {
		if (s < 0 || e >= ar.length) {
			throw new IndexOutOfBoundsException("Index out of Bounds, s=" + s + " > 0 and e=" + e + " < " + ar.length + ".");
		}
	}
	
	/**
	 * checks if the given {@code inds} are greater than 0
	 * @param inds
	 */
	public static void checkForGreaterZero(final int ... inds) {
		for (int i : inds) {
			if (i < 0) {
				throw new IllegalArgumentException("Indices must be greater than 0.");
			}
		}
	}
	
	/**
	 * checks if the given {@code inds} are greater than or equal to 0
	 * <br>(the input array is final)
	 * @param inds
	 */
	public static void checkForGreaterEqualZero(final int ... inds) {
		for (int i : inds) {
			if (i <= 0) {
				throw new IllegalArgumentException("Indices must be greater than or equal to0.");
			}
		}
	}
	
	public static void checkForGreaterEqualZero2(int ... inds) {
		for (int i : inds) {
			if (i <= 0) {
				throw new IllegalArgumentException("Indices must be greater than or equal to0.");
			}
		}
	}

	/**
	 * prints a 1D Object[] array
	 * @param ar
	 */
	public static void print(final Object[] ar) {
		System.out.println(Arrays.toString(ar));
	}
	
	/**
	 * prints a 1D double[] array
	 * @param ar
	 */
	public static void print(final double[] ar) {
		System.out.println(Arrays.toString(ar));
	}
	
	/**
	 * prints a 1D int[] array
	 * @param ar
	 */
	public static void print(final int[] ar) {
		System.out.println(Arrays.toString(ar));
	}
	
	/**
	 * prints the 2D object array
	 * @param ar
	 */
	public static void print(final Object[][] ar) {
		for (int i = 0; i < ar.length; i++) {
			System.out.println(Arrays.toString(ar[i]));
		}
	}
	
	/**
	 * prints the 2D double array
	 * <br>(Format corresponds with MATLAB matrix print out)
	 * @param matrix
	 */
	public static void print(final double[][] matrix) {
		System.out.println(toString(matrix));
	}
	
	/**
	 * returns the matrix as string
	 * <br>(Format corresponds with MATLAB matrix print out)
	 * @param matrix
	 * @return
	 */
	public static String toString(double[][] matrix) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < matrix.length; i++) {
			if (i == matrix.length - 1) {
				sb.append(" " + Arrays.toString(matrix[i]).replace("[", ""));
			} else if (i == 0) {
				sb.append(Arrays.toString(matrix[i]).replace("]", ";\n"));
			} else {
				sb.append(" " + Arrays.toString(matrix[i]).replace("[", "").replace("]", ";\n"));
			}
		}
		return sb.toString();
	}
	
	/**
	 * prints the 2D int array
	 * <br>(Format corresponds with MATLAB matrix print out)
	 * @param ar
	 */
	public static void print(final int[][] ar) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < ar.length; i++) {
			if (i == ar.length - 1) {
				sb.append(" " + Arrays.toString(ar[i]).replace("[", ""));
			} else if (i == 0) {
				sb.append(Arrays.toString(ar[i]).replace("]", ";\n"));
			} else {
				sb.append(" " + Arrays.toString(ar[i]).replace("[", "").replace("]", ";\n"));
			}
		}
		System.out.println(sb.toString());
	}
		
	public static String matrixDimensionsToString(double[][] X) {
		StringBuilder sb = new StringBuilder();
		int n1 = X[0].length;
		int m1 = X.length;
		sb.append("[").append(m1).append("]").append("[").append(n1).append("]");
		return sb.toString();
	}
}
