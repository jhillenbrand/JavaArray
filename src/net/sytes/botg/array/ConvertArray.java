package net.sytes.botg.array;

import java.lang.reflect.Array;

/**
 * 
 */

/**
 * @author jonas
 *
 */
public class ConvertArray {
	
	public static byte[] concat(byte[] b1, byte[] b2) {
		int len1 = b1.length;
	    int len2 = b2.length;
	    byte[] bn = new byte[len1 + len2];
	    System.arraycopy(b1, 0, bn, 0, len1);
	    System.arraycopy(b2, 0, bn, len1, len2);	
		return bn;
	}
	
	public static <T> T[] concat(T[] a, T[] b) {
	    int aLen = a.length;
	    int bLen = b.length;

	    @SuppressWarnings("unchecked")
	    T[] c = (T[]) Array.newInstance(a.getClass().getComponentType(), aLen + bLen);
	    System.arraycopy(a, 0, c, 0, aLen);
	    System.arraycopy(b, 0, c, aLen, bLen);

	    return c;
	}
	
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
