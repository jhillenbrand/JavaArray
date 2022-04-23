package net.sytes.botg.array.math;
/**
 * 
 */

/**
 * @author jonas
 *
 */
public class Round {
	
	public static double[] roundToDecimals(double[] x, int decimals) {
		double[] newX = new double[x.length];
		double scale = 0;
		int scaledXInt = 0;
		for (int i = 0; i < x.length; i++) {
			scale = Math.pow(10, decimals);
			scaledXInt = (int) x[i] * (int) scale;
			newX[i] = scaledXInt / scale;
		}
		return newX;
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
	
	public static double[] roundToDecimals(Number[] x, int decimals) {
		double[] newX = new double[x.length];
		double scale = 0;
		int scaledXInt = 0;
		for (int i = 0; i < x.length; i++) {
			scale = Math.pow(10, decimals);
			scaledXInt = (int) x[i] * (int) scale;
			newX[i] = scaledXInt / scale;
		}
		return newX;
	}
}
