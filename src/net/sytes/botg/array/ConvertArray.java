package net.sytes.botg.array;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.stream.LongStream;

/**
 * the unwrapping code is taken from org.apache.commons.lang3.ArrayUtils
 * https://www.apache.org/licenses/LICENSE-2.0
 * 
 * @author hillenbrand
 */
public class ConvertArray {
	
	/**
	 * concatenate two byte[] arrays
	 * @param b1
	 * @param b2
	 * @return
	 */
	public static byte[] concat(byte[] b1, byte[] b2) {
		int len1 = b1.length;
	    int len2 = b2.length;
	    byte[] bn = new byte[len1 + len2];
	    System.arraycopy(b1, 0, bn, 0, len1);
	    System.arraycopy(b2, 0, bn, len1, len2);	
		return bn;
	}
	
	/**
	 * concatenate two double[] arrays
	 * @param d1
	 * @param d2
	 * @return
	 */
	public static double[] concat(double[] d1, double[] d2) {
		int len1 = d1.length;
	    int len2 = d2.length;
	    double[] dn = new double[len1 + len2];
	    System.arraycopy(d1, 0, dn, 0, len1);
	    System.arraycopy(d2, 0, dn, len1, len2);	
		return dn;
	}
	
	/**
	 * concatenates two arrays
	 * @param <T>
	 * @param a
	 * @param b
	 * @return
	 */
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
	
	/**
	 * convert object[] to double[]
	 * @param objAr
	 * @return
	 */
	public static double[] toDouble(Object[] objAr) {
		double[] doubleAr = new double[objAr.length];
		for (int d = 0; d < doubleAr.length; d++) {
			doubleAr[d] = (double) objAr[d];
		}
		return doubleAr;
	}
	
	/**
	 * convert object[] to double[] using {@code Arrays.stream}
	 * @param objAr
	 * @return
	 */
	public static double[] toDouble2(Object[] objAr) {
		return Arrays.stream(objAr).mapToDouble(num -> (double) num).toArray();
	}
	
	/**
	 * convert int[] to double[]
	 * @param intAr
	 * @return
	 */
	public static double[] toDouble(int[] intAr) {
	    double[] doubleAr = new double[intAr.length];
	    for(int i=0; i < intAr.length; i++) {
	    	doubleAr[i] = intAr[i];
	    }
	    return doubleAr;
	}
	
	/**
	 * convert int[] to double[] using stream
	 * @param intAr
	 * @return
	 */
	public static double[] intToDouble2(int[] intAr) {
		return Arrays.stream(intAr).asDoubleStream().toArray();
	}
	
	public static double[] shortToDouble(short[] shortAr) {
		double[] doubleAr = new double[shortAr.length];
	    for(int i=0; i < shortAr.length; i++) {
	    	doubleAr[i] = shortAr[i];
	    }
	    return doubleAr;
	}
	
	/**
	 * convert long[] to double[]
	 * @param intAr
	 * @return
	 */
	public static double[] longToDouble(long[] longAr) {
	    double[] doubleAr = new double[longAr.length];
	    for(int i=0; i < longAr.length; i++) {
	    	doubleAr[i] = longAr[i];
	    }
	    return doubleAr;
	}
	
	/**
	 * convert long[] to double[] using stream
	 * @param longAr
	 * @return
	 */
	public static double[] longToDouble2(long[] longAr) {
		return Arrays.stream(longAr).asDoubleStream().toArray();
	}
	
	/**
     * <p>Converts an array of object Booleans to primitives.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Boolean} array, may be {@code null}
     * @return a {@code boolean} array, {@code null} if null array input
     */
    public static boolean[] unwrap(final Boolean[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final boolean[] result = new boolean[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].booleanValue();
        }
        return result;
    }

    /**
     * <p>Converts an array of object Booleans to primitives handling {@code null}.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Boolean} array, may be {@code null}
     * @param valueForNull  the value to insert if {@code null} found
     * @return a {@code boolean} array, {@code null} if null array input
     */
    public static boolean[] unwrap(final Boolean[] array, final boolean valueForNull) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final boolean[] result = new boolean[array.length];
        for (int i = 0; i < array.length; i++) {
            final Boolean b = array[i];
            result[i] = (b == null ? valueForNull : b.booleanValue());
        }
        return result;
    }

    /**
     * <p>Converts an array of object Bytes to primitives.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Byte} array, may be {@code null}
     * @return a {@code byte} array, {@code null} if null array input
     */
    public static byte[] unwrap(final Byte[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final byte[] result = new byte[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].byteValue();
        }
        return result;
    }

    /**
     * <p>Converts an array of object Bytes to primitives handling {@code null}.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Byte} array, may be {@code null}
     * @param valueForNull  the value to insert if {@code null} found
     * @return a {@code byte} array, {@code null} if null array input
     */
    public static byte[] unwrap(final Byte[] array, final byte valueForNull) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final byte[] result = new byte[array.length];
        for (int i = 0; i < array.length; i++) {
            final Byte b = array[i];
            result[i] = (b == null ? valueForNull : b.byteValue());
        }
        return result;
    }

    /**
     * <p>Converts an array of object Characters to primitives.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Character} array, may be {@code null}
     * @return a {@code char} array, {@code null} if null array input
     */
    public static char[] unwrap(final Character[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final char[] result = new char[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].charValue();
        }
        return result;
    }

    /**
     * <p>Converts an array of object Character to primitives handling {@code null}.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Character} array, may be {@code null}
     * @param valueForNull  the value to insert if {@code null} found
     * @return a {@code char} array, {@code null} if null array input
     */
    public static char[] unwrap(final Character[] array, final char valueForNull) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final char[] result = new char[array.length];
        for (int i = 0; i < array.length; i++) {
            final Character b = array[i];
            result[i] = (b == null ? valueForNull : b.charValue());
        }
        return result;
    }

   /**
     * <p>Converts an array of object Doubles to primitives.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Double} array, may be {@code null}
     * @return a {@code double} array, {@code null} if null array input
     */
    public static double[] unwrap(final Double[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final double[] result = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].doubleValue();
        }
        return result;
    }

    /**
     * <p>Converts an array of object Doubles to primitives handling {@code null}.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Double} array, may be {@code null}
     * @param valueForNull  the value to insert if {@code null} found
     * @return a {@code double} array, {@code null} if null array input
     */
    public static double[] unwrap(final Double[] array, final double valueForNull) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final double[] result = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            final Double b = array[i];
            result[i] = (b == null ? valueForNull : b.doubleValue());
        }
        return result;
    }
    
    public static Double[] wrap(final double[] array) {
    	if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final Double[] result = new Double[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i];
        }
        return result;
    }

    /**
     * <p>Converts an array of object Floats to primitives.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Float} array, may be {@code null}
     * @return a {@code float} array, {@code null} if null array input
     */
    public static float[] unwrap(final Float[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final float[] result = new float[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].floatValue();
        }
        return result;
    }

    /**
     * <p>Converts an array of object Floats to primitives handling {@code null}.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Float} array, may be {@code null}
     * @param valueForNull  the value to insert if {@code null} found
     * @return a {@code float} array, {@code null} if null array input
     */
    public static float[] unwrap(final Float[] array, final float valueForNull) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final float[] result = new float[array.length];
        for (int i = 0; i < array.length; i++) {
            final Float b = array[i];
            result[i] = (b == null ? valueForNull : b.floatValue());
        }
        return result;
    }

    /**
     * <p>Converts an array of object Integers to primitives.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Integer} array, may be {@code null}
     * @return an {@code int} array, {@code null} if null array input
     */
    public static int[] unwrap(final Integer[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final int[] result = new int[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].intValue();
        }
        return result;
    }

    /**
     * <p>Converts an array of object Integer to primitives handling {@code null}.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Integer} array, may be {@code null}
     * @param valueForNull  the value to insert if {@code null} found
     * @return an {@code int} array, {@code null} if null array input
     */
    public static int[] unwrap(final Integer[] array, final int valueForNull) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final int[] result = new int[array.length];
        for (int i = 0; i < array.length; i++) {
            final Integer b = array[i];
            result[i] = (b == null ? valueForNull : b.intValue());
        }
        return result;
    }

   /**
     * <p>Converts an array of object Longs to primitives.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Long} array, may be {@code null}
     * @return a {@code long} array, {@code null} if null array input
     */
    public static long[] unwrap(final Long[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final long[] result = new long[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].longValue();
        }
        return result;
    }

    /**
     * <p>Converts an array of object Long to primitives handling {@code null}.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Long} array, may be {@code null}
     * @param valueForNull  the value to insert if {@code null} found
     * @return a {@code long} array, {@code null} if null array input
     */
    public static long[] unwrap(final Long[] array, final long valueForNull) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final long[] result = new long[array.length];
        for (int i = 0; i < array.length; i++) {
            final Long b = array[i];
            result[i] = (b == null ? valueForNull : b.longValue());
        }
        return result;
    }

    /**
     * <p>Converts an array of object Shorts to primitives.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Short} array, may be {@code null}
     * @return a {@code byte} array, {@code null} if null array input
     */
    public static short[] unwrap(final Short[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final short[] result = new short[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].shortValue();
        }
        return result;
    }

    /**
     * <p>Converts an array of object Short to primitives handling {@code null}.
     *
     * <p>This method returns {@code null} for a {@code null} input array.
     *
     * @param array  a {@code Short} array, may be {@code null}
     * @param valueForNull  the value to insert if {@code null} found
     * @return a {@code byte} array, {@code null} if null array input
     */
    public static short[] unwrap(final Short[] array, final short valueForNull) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final short[] result = new short[array.length];
        for (int i = 0; i < array.length; i++) {
            final Short b = array[i];
            result[i] = (b == null ? valueForNull : b.shortValue());
        }
        return result;
    }	
}
