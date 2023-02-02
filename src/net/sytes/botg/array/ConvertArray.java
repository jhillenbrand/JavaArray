package net.sytes.botg.array;

import java.lang.reflect.Array;
import java.util.Arrays;

/**
 * the unwrapping code is taken from org.apache.commons.lang3.ArrayUtils
 * https://www.apache.org/licenses/LICENSE-2.0
 * 
 * @author hillenbrand
 */
public class ConvertArray {
	
	// Suppress default constructor for noninstantiability
	private ConvertArray() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
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
	public static double[] concat(double[] ar1, double[] ar2) {
		int len1 = ar1.length;
	    int len2 = ar2.length;
	    double[] ar = new double[len1 + len2];
	    System.arraycopy(ar1, 0, ar, 0, len1);
	    System.arraycopy(ar2, 0, ar, len1, len2);	
		return ar;
	}
	
	public static int[] concat(int[] ar1, int[] ar2) {
		int len1 = ar1.length;
	    int len2 = ar2.length;
	    int[] ar = new int[len1 + len2];
	    System.arraycopy(ar1, 0, ar, 0, len1);
	    System.arraycopy(ar2, 0, ar, len1, len2);	
		return ar;
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
	 * @param ar
	 * @return
	 */
	public static double[] toDouble(Object[] ar) {
		double[] doubleAr = new double[ar.length];
		for (int d = 0; d < doubleAr.length; d++) {
			doubleAr[d] = (double) ar[d];
		}
		return doubleAr;
	}
	
	/**
	 * converts the Double[] {@code ar} to Object[]
	 * @param ar
	 * @return
	 */
	public static Object[] toObjects(Double[] ar) {
		Object[] objAr = new Object[ar.length];
		for (int d = 0; d < objAr.length; d++) {
			objAr[d] = (Object) ar[d];
		}
		return objAr;
	}
	
	/**
	 * converts the double[] {@code ar} to Object[]
	 * @param ar
	 * @return
	 */
	public static Object[] toObjects(double[] ar) {
		Object[] objAr = new Object[ar.length];
		for (int d = 0; d < objAr.length; d++) {
			objAr[d] = (Object) ar[d];
		}
		return objAr;
	}
	
	/**
	 * converts an arbitrary object array into double array
	 * <br>by brute force parsing and casting the data types back to double
	 * @param objAr
	 * @return
	 */
	public static double[] parseToDouble(Object[] objAr) {
		double[] doubles = new double[objAr.length];
		int i = 0;
		for (Object o : objAr) {
			Class<?> clazz = o.getClass();
			if (clazz.equals(Double.class) || clazz.equals(double.class)) {
				doubles[i] = (double) o;
			} else if (clazz.equals(Integer.class) || clazz.equals(int.class)) {
				doubles[i] = (double) ((int) o);
			} else if (clazz.equals(Short.class) || clazz.equals(short.class)) {
				doubles[i] = (double) ((short) o);
			} else if (clazz.equals(Long.class) || clazz.equals(long.class) || clazz.equals(Long[].class) || clazz.equals(long[].class)) {
				doubles[i] = (double) ((long) o);
			} else if (clazz.equals(String.class)) {
				doubles[i] = Double.parseDouble(o.toString());
			} else if (clazz.equals(Float.class) || clazz.equals(float.class)) {
				doubles[i] = (double) ((float) o);
			} else if (clazz.equals(Boolean.class) || clazz.equals(boolean.class)) {
				doubles[i] = (((boolean) o) ? 1.0 : 0.0);
			} else {
				doubles[i] = Double.NaN;
			}
			++i;
		}
		return doubles;
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
	 * @param ar
	 * @return
	 */
	public static double[] longToDouble(long[] ar) {
	    double[] doubleAr = new double[ar.length];
	    for(int i=0; i < ar.length; i++) {
	    	doubleAr[i] = ar[i];
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
    
    /**
     * wraps an primitive double[] array into its Object wrapper
     * @param array
     * @return
     */
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
     * wraps an primitive int[] array into its Object wrapper
     * @param array
     * @return
     */
    public static Integer[] wrap(final int[] array) {
    	if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final Integer[] result = new Integer[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i];
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
    
    /**
     * unwraps {@code objects} into a 2D Object[][] array if elements are actually arrays
     * @param objects
     * @return
     */
    public static Object[][] unwrap(Object[] objects){
    	int m = objects.length;
    	Object o = objects[0];
    	if (o.getClass().isArray()) {
    		Object[] row0 = (Object[]) o;
    		int n = row0.length;
    		Object[][] matObjs = new Object[m][n];
    		for (int j = 0; j < m; j++) {
    			Object[] row = (Object[]) objects[j];
    			for (int k = 0; k < n; k++) {
    				matObjs[j][k] = row[k];
    			}
    		}
    		return matObjs;
    	} else {
    		return null;
    	}		
    }
    
}
