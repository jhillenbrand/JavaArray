package net.sytes.botg.array.math;

public class Scalar2Scalar {

	// Suppress default constructor for noninstantiability
	private Scalar2Scalar() {
		throw new AssertionError();
	}
		
	public static double roundToDecimals(double x, int decimals) {
		double scale = Math.pow(10, decimals);
		double scaledX = x * scale;
		int scaledXInt = (int) scaledX;
		double newX = scaledXInt / scale;
		return newX;
	}
	
	public static double roundToDecimals(Number x, int decimals) {
		double scale = Math.pow(10, decimals);
		double scaledX = x.doubleValue() * scale;
		long scaledXInt = (long) scaledX;
		double newX = scaledXInt / scale;
		return newX;
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
}
