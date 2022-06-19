package net.sytes.botg.array.math;

public class Vec2Mat {

	// Suppress default constructor for noninstantiability
	private Vec2Mat() {
		throw new AssertionError();
	}
	
	/**
	 * separates the given double array into windows of size windowSize
	 * <br>if dropUneven is true, then the remaining last ar.length % windowSize elements are dropped
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
	
}
