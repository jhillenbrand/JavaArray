package net.sytes.botg.array;

import java.util.Arrays;
import java.util.Random;

import net.sytes.botg.array.math.Scalar2Scalar;

public class ArrayUtility {

	// Suppress default constructor for noninstantiability
	private ArrayUtility() {
		throw new AssertionError();
	}
	
	private static final long BETTER_OF_AS_STREAM_SIZE = 100_000_000;
	
	public static void print2DArray(Object[][] ar) {
		for (int i = 0; i < ar.length; i++) {
			System.out.println(Arrays.toString(ar[i]));
		}
	}
	
	public static void print2DArray(double[][] ar) {
		for (int i = 0; i < ar.length; i++) {
			System.out.println(Arrays.toString(ar[i]));
		}
	}
	
	public static void print2DArray(int[][] ar) {
		for (int i = 0; i < ar.length; i++) {
			System.out.println(Arrays.toString(ar[i]));
		}
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
	
	public static double[] createRandomDoubleArray(long size) {
		double[] data;
		if (size < BETTER_OF_AS_STREAM_SIZE && size < Integer.MAX_VALUE) {
			data = new double[(int) size];
			Random r = new Random();
			for (int i = 0; i < size; i++) {
				data[i] = r.nextDouble();
			}
		} else {
			data = new Random().doubles(size).toArray();
		}
		return data;
	}
	
	public static int[] createRandomIntArray(long size) {
		int[] data;
		if (size < BETTER_OF_AS_STREAM_SIZE && size < Integer.MAX_VALUE) {
			data = new int[(int) size];
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
		double[] data = new double[size];
		for (int i = 0; i < size; i++) {
			data[i] = start + step * i;
		}
		return data;
	}
	
	/**
	 * returns a double array starting from {@code start} to {@code end} with equally spaced steps, so that the {@code size} is met
	 * @param start
	 * @param end
	 * @param size
	 * @return
	 */
	public static double[] linspace(double start, double end, int size) {
		double[] data = new double[size];
		double step = (end - start) / (size - 1);
		for (int i = 0; i < size; i++) {
			data[i] = start + step * i;
		}
		return data;
	}
	
	/**
	 * returns a double array starting from 0 and incrementing by 1 for {@code size} times
	 * @param size
	 * @return
	 */
	public static double[] linspace(int size) {
		return linspace(0, size, size);
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
		if (s >= 0 && e < ar.length) {
			double[] ar2 = new double[e - s];
			System.arraycopy(ar, s, ar2, 0, ar2.length);
			return ar2;
		} else {
			throw new IndexOutOfBoundsException("Index out of Range, s=" + s + " > 0 and e=" + e + " < " + ar.length + ".");
		}
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
		int m = (int) Scalar2Scalar.log2(n);		
		return m;
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
	
	public static void checkForEmpty(double[] ar) {
		if (ar.length == 0) {
			throw new IllegalArgumentException("ar must not be empty");
		}
	}
}
