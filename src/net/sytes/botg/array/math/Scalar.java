package net.sytes.botg.array.math;

import java.util.Random;

public class Scalar {

	// Suppress default constructor for noninstantiability
	private Scalar() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * --------------------------------------------------------------------------
	 * Scalar to Scalar
	 * --------------------------------------------------------------------------
	 */
	
	/**
	 * computes the norm (length) of the vector defined by [x, y]
	 * @param x
	 * @param y
	 * @return
	 */
	public static double norm(double x, double y) {
		return norm(x, y, 0.0);
	}
	
	/**
	 * computes the norm (length) of the vector defined by [x, y, z]
	 * @param x
	 * @param y
	 * @param z
	 * @return
	 */
	public static double norm(double x, double y, double z) {
		return Math.sqrt(x * x + y * y + z * z);
	}
	
	/**
	 * rounds a double to {@code decimals}
	 * @param x
	 * @param decimals
	 * @return
	 */
	public static double roundToDecimals(double x, int decimals) {
		double scale = Math.pow(10, decimals);
		double scaledX = x * scale;
		int scaledXInt = (int) scaledX;
		double newX = scaledXInt / scale;
		return newX;
	}
		
	/**
	 * rounds Number to {@code decimals} and returns double
	 * @param x
	 * @param decimals
	 * @return
	 */
	public static double roundToDecimals(Number x, int decimals) {
		double scale = Math.pow(10, decimals);
		double scaledX = x.doubleValue() * scale;
		long scaledXInt = (long) scaledX;
		double newX = scaledXInt / scale;
		return newX;
	}
	
	/**
	 * rounds the {@code x} to the next closest int
	 * <br>Example:
	 * <br>round(1.4) -> 1
	 * <br>round(1.5) -> 2
	 * @param x
	 * @return
	 */
	public static int round(double x) {
		int i = (int) x;
		double d = x - i;
		if (d >= 0.5) {
			return i + 1;
		} else {
			return i;
		}
	}
	
	/**
	 * computes the log of base 2 of @code d
	 * @param d
	 * @return
	 */
	public static double log2(double d) {
		return Math.log(d) / Math.log(2);
	}

	/**
	 * returns the signum for {@code d}
	 * @param d
	 * @return
	 */
	public static double sign(double d) {
		if (d < 0) {
			return -1;
		} else if (d > 0) {
			return 1;
		} else {
			return 0;
		}
	}
	
	/**
	 * updates existing mean value {@code mean} computed based on {@code n} samples with new datapoint {@code x}
	 * @return
	 */
	public static double updateMean(double mean, int n, double x) {
		return 1 / (n + 1) * (n * mean + x);
	}
	
	/**
	 * updates existing variance value {@code var} based on samples {@code n}, mean value {@code mean} and new data point {@code x}
	 * @param var
	 * @param n
	 * @param mean
	 * @param x
	 * @return
	 */
	public static double updateVariance(double var, int n, double mean, double x) {
		return 1 / (n + 1) * (n * var + Math.pow(x - mean, 2));
	}
	
	/**
	 * computes the next closest exponent for base 2 for {@code n}, that is smaller or equal to {@code n}
	 * @param n
	 * @return
	 */
	public static int nextLowerExponentForBase2(int n) {
		int m = (int) Scalar.log2(n);		
		return m;
	}
	
	/**
	 * computes the next greater exponent for base 2 for {@code n}, that is greater or equal to {@code n}
	 * @param n
	 * @return
	 */
	public static int nextGreaterExponentForBase2(int n) {
		int m = (int) Scalar.log2(n);
		if (n == Math.pow(2, m)) {
			return m;
		} else {
			return m + 1;
		}
	}
	
	/**
	 * --------------------------------------------------------------------------
	 * Scalar Generation Methods
	 * --------------------------------------------------------------------------
	 */
	
	/**
	 * returns a random int between {@code min} and {@code max}
	 * @param min
	 * @param max
	 * @return
	 */
	public static int randInt(int min, int max) {
		Scalar.checkForFirstSmallerSecond(min, max);
		Random r = new Random();
		return r.nextInt((max - min) + 1) + min;
	}

	/**
	 * create a random double between{@code min} and {@code max}
	 * @param min
	 * @param max
	 * @return
	 */
	public static double rand(double min, double max) {
		Random r = new Random();
		double randomValue = min + (max - min) * r.nextDouble();
		return randomValue;
	}
	
	/**
	 * --------------------------------------------------------------------------
	 * Scalar Utility Methods
	 * --------------------------------------------------------------------------
	 */
	
	/**
	 * If a number is divisible by 2 then it has its least significant bit (LSB) set to 0,
	 * 	if divisible by 4 then two LSB’s set to 0, if by 8 then three LSB’s set to 0 and so on.
	 * 	Keeping this in mind, a number n is divisible by 2m if (n &amp; ((1 &lt;&lt; m) – 1)) is equal to 0 else not.
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
	 * checks if the first argument is smaller than the second and throws {@code IllegalArgumentException} if true
	 * @param d1
	 * @param d2
	 */
	public static void checkForFirstSmallerSecond(double d1, double d2) {
		if (d1 >= d2) {
			throw new IllegalArgumentException("d1[" + d1 + "] must be smaller than d2[" + d2 + "]");
		}
	}
}
