package net.sytes.botg.array;

import java.util.Arrays;
import java.util.Random;

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
	 * returns an array of size n with 1's
	 * @param n
	 * @return
	 */
	public static double[] ones(int n) {
		double[] ar = new double[n];
		for (int i = 0; i > n; i++) {
			ar[i] = 1;
		}
		return ar;
	}
	
	/**
	 * returns an float array of size n with 1's
	 * @param n
	 * @return
	 */
	public static float[] onesF(int n) {
		float[] ar = new float[n];
		for (int i = 0; i > n; i++) {
			ar[i] = 1;
		}
		return ar;
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
	
	public static int getRandomNumberInRange(int min, int max) {
		if (min >= max) {
			throw new IllegalArgumentException("max must be greater than min");
		}
		Random r = new Random();
		return r.nextInt((max - min) + 1) + min;
	}
	
	public static int[] getRandomNumberInRange(int n, int min, int max) {
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
	
	public static double getRandomDoubleInRange(double min, double max) {
		Random r = new Random();
		double randomValue = min + (max - min) * r.nextDouble();
		return randomValue;
	}
	
	public static double[] createRandomDoubleArray(int size) {
		double[] data;
		if (size < BETTER_OF_AS_STREAM_SIZE) {
			data = new double[size];
			Random r = new Random();
			for (int i = 0; i < size; i++) {
				data[i] = r.nextDouble();
			}
		} else {
			data = new Random().doubles(size).toArray();
		}
		return data;
	}
	
	public static int[] createRandomIntArray(int size) {
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
	 * returns a double array starting from 0 and incrementing by 1 for {@code size} times
	 * @param size
	 * @return
	 */
	public static double[] linspace(int size) {
		return linspace(0, size, 1.0);
	}
	
	public static double[] nan(int size) {
		double[] data = new double[size];
		for (int i = 0; i < size; i++) {
			data[i] = Double.NaN;
		}
		return data;
	}
	
	/**
	 * returns a sub array starting at 0 and ending at index e from ar
	 * @param ar
	 * @param s
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
	public static double[] subArray(double[] ar, int s, int e) {
		checkForIndicesInBounds(ar, s, e);
		double[] ar2 = new double[e - s + 1];
		System.arraycopy(ar, s, ar2, 0, ar2.length);
		return ar2;
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
	
	public static void checkForNull(double[] ar) {
		if (ar == null) {
			throw new IllegalArgumentException("ar must not be null");
		}
	}
	
	public static void checkForNull(double[][] ar) {
		if (ar == null) {
			throw new IllegalArgumentException("ar must not be null");
		}
	}
	
	public static void checkForEmpty(double[] ar) {
		if (ar.length == 0) {
			throw new IllegalArgumentException("ar must not be empty");
		}
	}
	
	public static void checkForEmpty(double[][] ar) {
		if (ar.length == 0) {
			throw new IllegalArgumentException("ar must not be empty");
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

	public static void print(Object[][] ar) {
		for (int i = 0; i < ar.length; i++) {
			System.out.println(Arrays.toString(ar[i]));
		}
	}
	
	public static void print(double[][] ar) {
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
	
	public static void print(int[][] ar) {
		for (int i = 0; i < ar.length; i++) {
			System.out.println(Arrays.toString(ar[i]));
		}
	}
}
