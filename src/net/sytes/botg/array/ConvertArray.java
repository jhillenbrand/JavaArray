package net.sytes.botg.array;

/**
 * 
 */

/**
 * @author jonas
 *
 */
public class ConvertArray {
	
	/**
	 * This method converts an object array to a string. The elements are divided by a delimiter. 
	 * @param a
	 * @return
	 */
	public static String toString(Object[] a, String delimiter) {
        if (a == null) {
            return null;
        }
        int iMax = a.length - 1;
        if (iMax == -1)
            return null;

        StringBuilder b = new StringBuilder();
        for (int i = 0; ; i++) {
            b.append(String.valueOf(a[i]));
            if (i == iMax) {
                return b.toString();
            }
            b.append(delimiter);
        }
    }
	
	public static double[] copyFromIntToDoubleArray(int[] intAr) {
	    double[] doubleAr = new double[intAr.length];
	    for(int i=0; i < intAr.length; i++) {
	    	doubleAr[i] = intAr[i];
	    }
	    return doubleAr;
	}
	
	public static double[] copyFromLongToDoubleArray(long[] intAr) {
	    double[] doubleAr = new double[intAr.length];
	    for(int i=0; i < intAr.length; i++) {
	    	doubleAr[i] = intAr[i];
	    }
	    return doubleAr;
	}
}
