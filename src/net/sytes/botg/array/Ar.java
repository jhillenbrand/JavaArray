package net.sytes.botg.array;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import net.sytes.botg.array.math.Vec;

/**
 * the unwrapping code is taken from org.apache.commons.lang3.ArrayUtils
 * https://www.apache.org/licenses/LICENSE-2.0
 * 
 * @author hillenbrand
 */
public class Ar {
	
	// Suppress default constructor for noninstantiability
	private Ar() {
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
	 * appends a {@code prefix} inf front of all elements of {@code stringAr}
	 * @param stringAr
	 * @return
	 */
	public static String[] prefix(String[] stringAr, String prefix) {
		String[] newAr = new String[stringAr.length];
		for (int s = 0; s < stringAr.length; s++) {
			newAr[s] = prefix + stringAr[s];
		}
		return newAr;
	}
	
	/**
	 * appends a {@code postfix} to all elements of {@code stringAr}
	 * @param stringAr
	 * @return
	 */
	public static String[] postfix(String[] stringAr, String postfix) {
		String[] newAr = new String[stringAr.length];
		for (int s = 0; s < stringAr.length; s++) {
			newAr[s] = stringAr[s] + postfix;
		}
		return newAr;
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
     * wraps an primitive float[] array into its Object wrapper
     * @param array
     * @return
     */
    public static Float[] wrap(final float[] array) {
    	if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final Float[] result = new Float[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i];
        }
        return result;
    }
    
    /**
     * wraps an primitive double[] array into its Object wrapper
     * @param array
     * @return
     */
    public static Long[] wrap(final long[] array) {
    	if (array == null) {
            return null;
        } else if (array.length == 0) {
            return null;
        }
        final Long[] result = new Long[array.length];
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
    
	/**
	 * delete the element of double array {@code x} at index {@code i} (reducing the array length by one), by returning a copy of original array except this index {@code i}
	 * @param i
	 * @return
	 */
	public static double[] removeAt(double[] x, int i) {
		double[] y = new double[x.length - 1];
		int j = 0;
		for (int k = 0; k < x.length; k++) {
			if (k != i) {
				y[j] = x[k];
				++j;
			}
		}
		return y;
	}
	    
	/**
	 * prints a 1D Object[] array
	 * @param ar
	 */
	public static void print(final Object[] ar) {
		System.out.println(Arrays.toString(ar));
	}
	
	/**
	 * prints a 1D double[] array
	 * @param ar
	 */
	public static void print(final double[] ar) {
		System.out.println(Arrays.toString(ar));
	}
	
	/**
	 * prints a 1D int[] array
	 * @param ar
	 */
	public static void print(final int[] ar) {
		System.out.println(Arrays.toString(ar));
	}
		
	/**
	 * copies the specified array {@code x}
	 * @param x
	 * @return
	 */
	public static double[] copy(double[] x) {
		return x.clone();
	}
	
	/**
	 * copies the specified array {@code x}
	 * @param x
	 * @return
	 */
	public static byte[] copy(byte[] x) {
		return x.clone();
	}
	
	/**
	 * copies the specified array {@code x}
	 * @param x
	 * @return
	 */
	public static int[] copy(int[] x) {
		return x.clone();
	}
	
	/**
	 * copies the specified array {@code x}
	 * @param x
	 * @return
	 */
	public static long[] copy(long[] x) {
		return x.clone();
	}
	
	/**
	 * copies the specified array {@code x}
	 * @param x
	 * @return
	 */
	public static float[] copy(float[] x) {
		return x.clone();
	}
	
	/**
	 * copies the specified array {@code x}
	 * @param x
	 * @return
	 */
	public static String[] copy(String[] x) {
		return x.clone();
	}
	
	/**
	 * copies the specified array {@code x}
	 * @param x
	 * @return
	 */
	public static Object[] copy(Object[] x) {
		return x.clone();
	}
	
	/**
	 * returns true if set specified object {@code o} is an array
	 * @param o
	 * @return
	 */
	public static boolean isArray(Object o) {
		if (o != null && o.getClass().isArray()) {
			return true;
		} else {
			return false;
		}
	}
	
	public static void checkForEmpty(double[] ar) {
		if (ar.length == 0) {
			throw new IllegalArgumentException("ar must not be empty");
		}
	}
	
	public static void checkForEmpty(int[] ar) {
		if (ar.length == 0) {
			throw new IllegalArgumentException("ar must not be empty");
		}
	}
	
	/**
	 * checks whether array is empty
	 * @param ar
	 */
	public static void checkForEmpty(float[] ar) {
		if (ar.length == 0) {
			throw new IllegalArgumentException("ar must not be empty");
		}
	}
	
	/**
	 * checks whether matrix is empty
	 * @param matrix
	 */
	public static void checkForEmpty(float[][] matrix) {
		if (matrix.length == 0) {
			throw new IllegalArgumentException("matrix must not be empty");
		}
	}
	
	public static void checkForEmpty2(double[] ... x) {
		for (int i = 0; i < x.length; i++) {
			if (x[i].length == 0) {
				throw new IllegalArgumentException("input arrays must not be empty");
			}
		}
	}	
		
	/**
	 * checks if the given {@code inds} are greater than 0
	 * @param inds
	 */
	public static void checkForGreaterZero(final int ... inds) {
		for (int i : inds) {
			if (i < 0) {
				throw new IllegalArgumentException("Indices must be greater than 0.");
			}
		}
	}
	
	/**
	 * checks if the given {@code inds} are greater than or equal to 0
	 * <br>(the input array is final)
	 * @param inds
	 */
	public static void checkForGreaterEqualZero(final int[] inds) {
		for (int i : inds) {
			if (i <= 0) {
				throw new IllegalArgumentException("Indices must be greater than or equal to 0.");
			}
		}
	}
	
	public static void checkForGreaterEqualZero2(int[] inds) {
		for (int i : inds) {
			if (i <= 0) {
				throw new IllegalArgumentException("Indices must be greater than or equal to 0.");
			}
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
	
	/**
	 * checks if the start and end indices are within array bounds
	 * @param ar
	 * @param s
	 * @param e
	 */
	public static void checkForIndicesInBounds(int[] ar, int s, int e) {
		if (s < 0 || e >= ar.length) {
			throw new IndexOutOfBoundsException("Index out of Bounds, s=" + s + " > 0 and e=" + e + " < " + ar.length + ".");
		}
	}
	
	/**
	 * checks if the start and end indices are within array bounds
	 * @param ar
	 * @param s
	 * @param e
	 */
	public static void checkForIndicesInBounds(String[] ar, int s, int e) {
		if (s < 0 || e >= ar.length) {
			throw new IndexOutOfBoundsException("Index out of Bounds, s=" + s + " > 0 and e=" + e + " < " + ar.length + ".");
		}
	}
	
	/**
	 * check if argument not NULL
	 * @param ar
	 */
	public static void checkForNull(double[] ar) {
		if (ar == null) {
			throw new IllegalArgumentException("ar must not be null");
		}
	}
	
	/**
	 * check if argument not NULL
	 * @param ar
	 */
	public static void checkForNull(int[] ar) {
		if (ar == null) {
			throw new IllegalArgumentException("ar must not be null");
		}
	}
	
	/**
	 * check if argument not NULL
	 * @param ar
	 */
	public static void checkForNull(float[] ar) {
		if (ar == null) {
			throw new IllegalArgumentException("ar must not be null");
		}
	}
	
	/**
	 * check if arguments are not NULL
	 * @param ar
	 */
	public static void checkForNull2(double[] ... ar) {
		for (int a = 0; a < ar.length; a++) {
			if (ar[a] == null) {
				throw new IllegalArgumentException("arrays must not be NULL");
			}
		}
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
	 * checks whether specified {@code arrays} have same length
	 * @param arrays
	 */
	public static void checkForEqualDimensions(double[] ... arrays) {
		if (arrays.length == 0) {
			throw new IllegalArgumentException("specified arrays were empty");
		}		
		int n = arrays[0].length;
		for (double[] ar : arrays) {
			if (ar.length != n) {
				throw new IllegalArgumentException("not all specified arrays have the same length");
			}
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
	 * checks whether there are at least {@code n} elements in {@code ar}
	 * @param ar
	 * @param n
	 */
	public static void checkForAtLeastNElements(double[] ar, int n) {
		if (ar.length < n) {
			throw new IllegalArgumentException("ar must have at least " + n + " elements");
		}
	}
	
	/**
	 * check if both matrices have same dimensions
	 * @param X
	 * @param Y
	 */
	public static void checkForMatchingDimensions(double[][] X, double[][] Y) {
		if (X.length != Y.length || X[0].length != Y[0].length) {
			throw new IllegalArgumentException("Dimension mismatch of matrices");
		}
	}
	
	/**
	 * checks if the matrix passed is null
	 * @param X
	 * @throws IllegalArgumentException when the matrix passed is null
	 */
	public static void checkForNull(double[][] X) {
		if (X == null) {
			throw new IllegalArgumentException("matrix must not be null");
		}
	}
	
	/**
	 * checks if the matrix passed is null
	 * @param X
	 * @throws IllegalArgumentException when the matrix passed is null
	 */
	public static void checkForNull(float[][] X) {
		if (X == null) {
			throw new IllegalArgumentException("matrix must not be null");
		}
	}
	
	/**
	 * checks if the matrix passed is null
	 * @param X
	 * @throws IllegalArgumentException when the matrix passed is null
	 */
	public static void checkForNull(long[][] X) {
		if (X == null) {
			throw new IllegalArgumentException("matrix must not be null");
		}
	}
	
	/**
	 * checks if the matrix passed is null
	 * @param X
	 * @throws IllegalArgumentException when the matrix passed is null
	 */
	public static void checkForNull(int[][] X) {
		if (X == null) {
			throw new IllegalArgumentException("matrix must not be null");
		}
	}
	
	/**
	 * checks if the matrix passed is null
	 * @param X
	 * @throws IllegalArgumentException when the matrix passed is null
	 */
	public static void checkForNull(Object[][] X) {
		if (X == null) {
			throw new IllegalArgumentException("matrix must not be null");
		}
	}
	
	/**
	 * checks if the matrix passed is null
	 * @param X
	 * @throws IllegalArgumentException when the matrix passed is null
	 */
	public static void checkForNull(byte[][] X) {
		if (X == null) {
			throw new IllegalArgumentException("matrix must not be null");
		}
	}
	
	/**
	 * checks if none of the matrices passed are null
	 * @param X
	 * @throws IllegalArgumentException when one of the matrices passed is null
	 */
	public static void checkForNull(double[][] ... X) {
		for (int i = 0; i < X.length; i++) {
			if (X[i] == null) {
				throw new IllegalArgumentException("matrices must not be null");
			}
		}
	}
	
	public static void checkForEmpty(double[][] X) {
		if (X.length == 0) {
			throw new IllegalArgumentException("matrix X must not be empty");
		} else {
			for (int i = 0; i < X.length; i++) {
				if (X[i].length == 0) {
					throw new IllegalArgumentException("matrix X must not contain empty vectors");
				}
			}
		}
	}
	
	public static void checkForEmpty(long[][] X) {
		if (X.length == 0) {
			throw new IllegalArgumentException("matrix X must not be empty");
		} else {
			for (int i = 0; i < X.length; i++) {
				if (X[i].length == 0) {
					throw new IllegalArgumentException("matrix X must not contain empty vectors");
				}
			}
		}
	}
	
	public static void checkForEmpty(int[][] X) {
		if (X.length == 0) {
			throw new IllegalArgumentException("matrix X must not be empty");
		} else {
			for (int i = 0; i < X.length; i++) {
				if (X[i].length == 0) {
					throw new IllegalArgumentException("matrix X must not contain empty vectors");
				}
			}
		}
	}
	
	public static void checkForEmpty(double[][] ... X) {
		for (int i = 0; i < X.length; i++) {
			if (X[i].length == 0) {
				throw new IllegalArgumentException("ar must not be empty");
			}
		}
	}

	public static void checkForEmpty(Object[][] X) {
		if (X.length == 0) {
			throw new IllegalArgumentException("matrix X must not be empty");
		} else {
			for (int i = 0; i < X.length; i++) {
				if (X[i].length == 0) {
					throw new IllegalArgumentException("matrix X must not contain empty vectors");
				}
			}
		}
	}
	
	public static void checkForEmpty(byte[][] X) {
		if (X.length == 0) {
			throw new IllegalArgumentException("matrix X must not be empty");
		} else {
			for (int i = 0; i < X.length; i++) {
				if (X[i].length == 0) {
					throw new IllegalArgumentException("matrix X must not contain empty vectors");
				}
			}
		}
	}
	
	/**
	 * prints the 2D object array
	 * @param ar
	 */
	public static void print(final Object[][] ar) {
		for (int i = 0; i < ar.length; i++) {
			System.out.println(Arrays.toString(ar[i]));
		}
	}
	
	/**
	 * prints the 2D double array
	 * <br>(Format corresponds with MATLAB matrix print out)
	 * @param matrix
	 */
	public static void print(final double[][] matrix) {
		System.out.println(toString(matrix));
	}
	
	/**
	 * prints the 2D float array
	 * <br>(Format corresponds with MATLAB matrix print out)
	 * @param matrix
	 */
	public static void print(final float[][] matrix) {
		System.out.println(toString(matrix));
	}
	
	/**
	 * deep clone an matrix {@code X}
	 * @param X
	 * @return
	 */
	public static double[][] copy(double[][] X){
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		int r = X.length;
		double[][] Y = new double[r][];
		for (int i = 0; i < r; i++) {
			Y[i] = X[i].clone();
		}
		return Y;
	}
	
	/**
	 * deep clone an matrix {@code X} using loop
	 * @param X
	 * @return
	 */
	public static double[][] copy2(double[][] X){
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		int r = X.length;
		int c = X[0].length;
		double[][] Y = new double[r][c];
		for (int i = 0; i < c; i++) {
			for (int j = 0; j < r; j++) {
				Y[j][i] = X[j][i];
			}
		}
		return Y;
	}
	
	public static float[][] copy(float[][] X){
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		int r = X.length;
		float[][] Y = new float[r][];
		for (int i = 0; i < r; i++) {
			Y[i] = X[i].clone();
		}
		return Y;
	}
	
	public static int[][] copy(int[][] X){
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		int r = X.length;
		int[][] Y = new int[r][];
		for (int i = 0; i < r; i++) {
			Y[i] = X[i].clone();
		}
		return Y;
	}	
	
	public static long[][] copy(long[][] X){
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		int r = X.length;
		long[][] Y = new long[r][];
		for (int i = 0; i < r; i++) {
			Y[i] = X[i].clone();
		}
		return Y;
	}
	
	public static Object[][] copy(Object[][] X){
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		int r = X.length;
		Object[][] Y = new Object[r][];
		for (int i = 0; i < r; i++) {
			Y[i] = X[i].clone();
		}
		return Y;
	}
	
	public static byte[][] copy(byte[][] X){
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		int r = X.length;
		byte[][] Y = new byte[r][];
		for (int i = 0; i < r; i++) {
			Y[i] = X[i].clone();
		}
		return Y;
	}
		
	/**
	 * returns the matrix as string
	 * <br>(Format corresponds with MATLAB matrix print out)
	 * @param matrix
	 * @return
	 */
	public static String toString(double[][] matrix) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < matrix.length; i++) {
			if (i == matrix.length - 1 && matrix.length != 1) {
				sb.append(" " + Arrays.toString(matrix[i]).replace("[", ""));
			} else if (matrix.length == 1) {
				sb.append(Arrays.toString(matrix[i]));
			} else if (i == 0) {
				sb.append(Arrays.toString(matrix[i]).replace("]", ";\n"));
			} else {
				sb.append(" " + Arrays.toString(matrix[i]).replace("[", "").replace("]", ";\n"));
			}
		}
		return sb.toString();
	}
	
	/**
	 * returns the matrix as string
	 * <br>(Format corresponds with MATLAB matrix print out)
	 * @param matrix
	 * @return
	 */
	public static String toString(float[][] matrix) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < matrix.length; i++) {
			if (i == matrix.length - 1 && matrix.length != 1) {
				sb.append(" " + Arrays.toString(matrix[i]).replace("[", ""));
			} else if (matrix.length == 1) {
				sb.append(Arrays.toString(matrix[i]));
			} else if (i == 0) {
				sb.append(Arrays.toString(matrix[i]).replace("]", ";\n"));
			} else {
				sb.append(" " + Arrays.toString(matrix[i]).replace("[", "").replace("]", ";\n"));
			}
		}
		return sb.toString();
	}
	
	/**
	 * returns the matrix as string
	 * <br>(Format corresponds with MATLAB matrix print out)
	 * @param matrix
	 * @return
	 */
	public static String toString(long[][] matrix) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < matrix.length; i++) {
			if (i == matrix.length - 1 && matrix.length != 1) {
				sb.append(" " + Arrays.toString(matrix[i]).replace("[", ""));
			} else if (matrix.length == 1) {
				sb.append(Arrays.toString(matrix[i]));
			} else if (i == 0) {
				sb.append(Arrays.toString(matrix[i]).replace("]", ";\n"));
			} else {
				sb.append(" " + Arrays.toString(matrix[i]).replace("[", "").replace("]", ";\n"));
			}
		}
		return sb.toString();
	}
	
	/**
	 * returns the matrix as string
	 * <br>(Format corresponds with MATLAB matrix print out)
	 * @param matrix
	 * @return
	 */
	public static String toString(int[][] matrix) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < matrix.length; i++) {
			if (i == matrix.length - 1 && matrix.length != 1) {
				sb.append(" " + Arrays.toString(matrix[i]).replace("[", ""));
			} else if (matrix.length == 1) {
				sb.append(Arrays.toString(matrix[i]));
			} else if (i == 0) {
				sb.append(Arrays.toString(matrix[i]).replace("]", ";\n"));
			} else {
				sb.append(" " + Arrays.toString(matrix[i]).replace("[", "").replace("]", ";\n"));
			}
		}
		return sb.toString();
	}
	
	/**
	 * prints the 2D int array
	 * <br>(Format corresponds with MATLAB matrix print out)
	 * @param ar
	 */
	public static void print(final int[][] ar) {
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
	
	/**
	 * method returns true/false whether any of the array elements ar are contained in String s
	 * @param s
	 * @param ar
	 * @return 
	 */
	public static boolean contains(double[] x, double d) {
		for (double xi : x) {
			if (xi == d) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * method returns true/false whether any of the array elements ar are contained in String s
	 * @param s
	 * @param ar
	 * @return 
	 */
	public static boolean contains(int[] x, int i) {
		for (int xi : x) {
			if (xi == i) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * method returns true/false whether any of the array elements {@code x} contain String {@code s}
	 * @param s
	 * @param ar
	 * @return 
	 */
	public static boolean contains(String[] x, String s) {
		for (String xi : x) {
			if (xi.contains(s)) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * method returns true/false whether any of the array elements {@code ar} contain Object {@code e}
	 * the method is typed during runtime
	 * @param <T>
	 * @param ar
	 * @param e
	 * @return
	 */
	public static <T> boolean contains(T[] ar, T e) {
		for (T a : ar) {
			if (a.equals(e)) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * counts the occurences of {@code s} in array {@code ar}
	 * if {@code caseSensitive} is specified the search neglects capital or lowercase letters
	 * @param ar
	 * @param s
	 * @param caseSensitive
	 * @return
	 */
    public static int countOccurences(String[] ar, String s, boolean caseSensitive) {
    	int count = -1;
        for (int i = 0; i < ar.length; i++) {
        	if (!caseSensitive) {
	        	if (s.toLowerCase().equals(ar[i].toLowerCase())){
	        		count++;
	            }
        	} else {
        		if (s.equals(ar[i])){
	        		count++;
	            }
        	}
        }
        return count;
    }   
    	
	/**
	 * returns a sub array starting at index s and ending at index e from ar
	 * @param ar
	 * @param s
	 * @param e
	 * @return
	 */
	public static double[] sub(final double[] ar, int s, int e) {
		Ar.checkForIndicesInBounds(ar, s, e);
		double[] ar2 = new double[e - s + 1];
		System.arraycopy(ar, s, ar2, 0, ar2.length);
		return ar2;
	}
	
	/**
	 * returns a sub array starting at 0 and ending at index e from ar
	 * @param ar
	 * @param e
	 * @return
	 */
	public static double[] sub(double[] ar, int e) {
		return sub(ar, 0, e);
	}
	
	/**
	 * returns a sub array of {@code x} between {@code start} and {@code end} with {@code step}
	 * @param x
	 * @param start
	 * @param end
	 * @param step
	 * @return
	 */
	public static double[] sub(double[] x, int start, int end, int step) {
		int n = (end - start + 1) / step;
		double[] y = new double[n];
		int j = 0;
		for (int i = start; i < end; i = i + step) {
			y[j] = x[i];
			++j;
		}
		return y;
	}
	
	/**
	 * returns a sub array starting at 0 and ending at index e from ar
	 * @param ar
	 * @param e
	 * @return
	 */
	public static int[] sub(int[] ar, int e) {
		return sub(ar, 0, e);
	}
	
	/**
	 * returns a sub array starting at index s and ending at index e from ar
	 * @param ar
	 * @param s
	 * @param e
	 * @return
	 */
	public static int[] sub(final int[] ar, int s, int e) {
		Ar.checkForIndicesInBounds(ar, s, e);
		int[] ar2 = new int[e - s + 1];
		System.arraycopy(ar, s, ar2, 0, ar2.length);
		return ar2;
	}
	
	/**
	 * returns a sub array starting at index s and ending at index e from ar
	 * @param ar
	 * @param s
	 * @param e
	 * @return
	 */
	public static String[] sub(final String[] ar, int s, int e) {
		Ar.checkForIndicesInBounds(ar, s, e);
		String[] ar2 = new String[e - s + 1];
		System.arraycopy(ar, s, ar2, 0, ar2.length);
		return ar2;
	}
	
	/**
	 * removes {@code start} samples from the beginning of the {@code ar} and {@code end} samples from the end of {@code ar} 
	 * <br>and returns the result as new double[] array 
	 * @param start
	 * @param end
	 * @return
	 */
	public static double[] trim(double[] ar, int start, int end) {
		Ar.checkForGreaterEqualZero(new int[] {start, end});
		int n = ar.length;
		if (start + end > n) {
			throw new IllegalArgumentException("sum of start and end must be smaller than length of ar");
		}
		return Ar.sub(ar, start, n - end); 
	}
		
	/**
	 * returns the elements specified with {@code inds} in {@code data} into new array
	 * @param inds
	 * @param data
	 * @return double[]
	 */
	public static double[] elementsAt(double[] data, boolean[] inds) {
		if (data.length != inds.length) {
			throw new IllegalArgumentException("length of data and inds array must have same length");
		}
		int n = Vec.sum(inds);
		double[] newData = new double[n];
		int j = 0;
		for (int i = 0; i < data.length; i++) {
			if (inds[i]) {
				newData[j] = data[i];
				++j;
			}
		}
		return newData;
	}
	
	/**
	 * returns the elements at indices {@code inds} in {@code ar} into new array
	 * @param inds
	 * @param data
	 * @return double[]
	 */
	public static double[] elementsAt(double[] ar, int[] inds) {
		int n = inds.length;
		double[] ar2 = new double[n];
		for (int i = 0; i < n; i++) {
			ar2[i] = ar[inds[i]];
		}
		return ar2;
	}
	
	/**
	 * returns the elements at indices {@code inds} in {@code ar} into new array
	 * @param inds
	 * @param data
	 * @return double[]
	 */
	public static double[] elementsAt(double[] ar, double[] inds) {
		int n = inds.length;
		double[] ar2 = new double[n];
		for (int i = 0; i < n; i++) {
			ar2[i] = ar[(int) inds[i]];
		}
		return ar2;
	}
	
	/**
	 * returns the elements at indices {@code inds} in {@code ar} into new array
	 * @param inds
	 * @param data
	 * @return int[]
	 */
	public static int[] elementsAt(int[] ar, int[] inds) {
		int n = inds.length;
		int[] ar2 = new int[n];
		for (int i = 0; i < n; i++) {
			ar2[i] = ar[inds[i]];
		}
		return ar2;
	}

	/**
	 * add the element x to the end of the array ar and return the new array
	 * @param <T> generic type
	 * @param x new sample
	 * @param ar existing array
	 * @return generic array T[] 
	 */
	public static <T> T[] append(T[] ar, T x) {
		T[] newAr = Arrays.copyOf(ar, ar.length + 1);
		newAr[ar.length] = x;
		return newAr;
	}
	
	/**
	 * add double d at the end of array and return new array
	 * @param d
	 * @param ar
	 * @return
	 */
	public static double[] append(double[] ar, double d) {
		double[] newAr = new double[ar.length + 1];
		System.arraycopy(ar, 0, newAr, 0, ar.length);
		newAr[ar.length] = d;
		return newAr;
	}
	
	/**
	 * returns an array with unique elements based on values in {@code x}
	 * @param x
	 * @return
	 */
	public static double[] unique(final double[] x) {
		Ar.checkForNull(x);
		Ar.checkForEmpty(x);
		int n = x.length;
		HashMap<Double, Double> xu = new LinkedHashMap<Double, Double>();
		//HashMap<Double, Double> xu = new HashMap<Double, Double>();
		for (int i = 0; i < n; i++) {
			xu.put(x[i], x[i]);
		}
		return Ar.unwrap(xu.values().toArray(new Double[xu.size()]));
	}
	
	/**
	 * returns an array with unique elements based on values of {@code x}
	 * in ascending numerical order
	 * @param x
	 * @return
	 */
	public static int[] unique(final int[] x) {
		Ar.checkForNull(x);
		Ar.checkForEmpty(x);
		int n = x.length;
		HashMap<Integer, Integer> xu = new LinkedHashMap<Integer, Integer>();
		//HashMap<Double, Double> xu = new HashMap<Double, Double>();
		for (int i = 0; i < n; i++) {
			xu.put(x[i], x[i]);
		}
		int[] index = Ar.unwrap(xu.values().toArray(new Integer[xu.size()]));
		Arrays.sort(index);
		return index;
	}
		
	/**
	 * returns an matrix with unique elements based on values in {@code X}
	 * @param X
	 * @return
	 */
	public static double[][] unique(final double[][] X) {
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		int n = X.length;
		HashMap<String, double[]> XU = new LinkedHashMap<String, double[]>();
		//HashMap<Double[], Double[]> XU = new HashMap<Double[], Double[]>();
		for (int i = 0; i < n; i++) {
			XU.put(Arrays.toString(X[i]), X[i]);
		}
		Object[] objs = XU.values().toArray();
		
		double[][] Y = new double[objs.length][];
		for (int i = 0; i < Y.length; i++) {
			Y[i] = (double[]) objs[i];
		}
		
		return Y;
	}

	/**
	 * searches the first index in {@code x}, where the element of {@code x} is equal to {@code d}
	 * <br>if no element is found -1 is returned
	 * @param x
	 * @param d
	 * @return
	 */
	public static int findFirst(double[] x, double d) {
		for (int i = 0; i < x.length; i++) {
			if (x[i] == d) {
				return i;
			}
		}
		return -1;
	}
	
	/**
	 * searches the indices in {@code x}, where the element of {@code x} is equal to {@code d}
	 * <br>if no element is found NULL is returned
	 * @param x
	 * @param d
	 * @return
	 */
	public static int[] find(double[] x, double d) {
		boolean none = true;
		int[] inds = new int[x.length];
		int c = 0;
		for (int i = 0; i < x.length; i++) {
			if (x[i] == d) {
				inds[c] = i;
				++c;
				none = false;
			}
		}
		if (none) {
			return null;
		}
		inds = Ar.sub(inds, c - 1);
		return inds;
	}
	
	/**
	 * searches the indices in {@code x}, where the element of {@code x} is equal to {@code xi}
	 * <br>if no element is found NULL is returned
	 * @param x
	 * @param xi
	 * @return
	 */
	public static int[] find(int[] x, int xi) {
		boolean none = true;
		int[] inds = new int[x.length];
		int c = 0;
		for (int i = 0; i < x.length; i++) {
			if (x[i] == xi) {
				inds[c] = i;
				++c;
				none = false;
			}
		}
		if (none) {
			return null;
		}
		inds = Ar.sub(inds, c - 1);
		return inds;
	}

	/**
	 * returns the index of {@code str} found in {@code ar}
	 * <br>if {@code ar} does not contain {@code str}, then -1 is returned 
	 * @param ar
	 * @param str
	 * @return
	 */
	public static int indexAt(final String[] ar, final String str) {
		int i = -1;
		for (int s = 0; s < ar.length; s++) {
			if (ar[s].contentEquals(str)) {
				i = s;
				break;
			}
		}
		return i;
	}
	
	/**
	 * checks whether all elements in {@code ar} are the same
	 * @param <T>
	 * @param ar
	 * @return
	 */
	public static <T> boolean areAllEqual(T[] ar) {
		if (ar.length > 1) {
			T elem = ar[0];
			for (int i = 1; i < ar.length; i++) {
				if (!ar[i].equals(elem)) {
					return false;
				}
				elem = ar[i];
			}
			return true;
		} else {
			return true;
		}
	}
	
	/**
	 * checks whether all elements in {@code ar} are the same
	 * @param ar
	 * @return
	 */
	public static boolean areAllEqual(int[] ar) {
		if (ar.length > 1) {
			int elem = ar[0];
			for (int i = 1; i < ar.length; i++) {
				if (ar[i] != elem) {
					return false;
				}
				elem = ar[i];
			}
			return true;
		} else {
			return true;
		}
	}
	
	/**
	 * checks whether all elements in {@code ar} are the same
	 * @param ar
	 * @return
	 */
	public static boolean areAllEqual(double[] ar) {
		if (ar.length > 1) {
			double elem = ar[0];
			for (int i = 1; i < ar.length; i++) {
				if (ar[i] != elem) {
					return false;
				}
				elem = ar[i];
			}
			return true;
		} else {
			return true;
		}
	}
	
	/**
	 * checks whether all elements in {@code ar} are the same
	 * @param ar
	 * @return
	 */
	public static boolean areAllEqual(float[] ar) {
		if (ar.length > 1) {
			float elem = ar[0];
			for (int i = 1; i < ar.length; i++) {
				if (ar[i] != elem) {
					return false;
				}
				elem = ar[i];
			}
			return true;
		} else {
			return true;
		}
	}
	
	/**
	 * checks whether all elements in {@code ar} are the same
	 * @param ar
	 * @return
	 */
	public static boolean areAllEqual(long[] ar) {
		if (ar.length > 1) {
			long elem = ar[0];
			for (int i = 1; i < ar.length; i++) {
				if (ar[i] != elem) {
					return false;
				}
				elem = ar[i];
			}
			return true;
		} else {
			return true;
		}
	}
	
	/**
	 * checks whether all elements in {@code ar} are the same
	 * @param ar
	 * @return
	 */
	public static boolean areAllEqual(byte[] ar) {
		if (ar.length > 1) {
			byte elem = ar[0];
			for (int i = 1; i < ar.length; i++) {
				if (ar[i] != elem) {
					return false;
				}
				elem = ar[i];
			}
			return true;
		} else {
			return true;
		}
	}
	
	/**
	 * checks whether all elements in {@code ar} are the same
	 * @param ar
	 * @return
	 */
	public static boolean areAllEqual(short[] ar) {
		if (ar.length > 1) {
			short elem = ar[0];
			for (int i = 1; i < ar.length; i++) {
				if (ar[i] != elem) {
					return false;
				}
				elem = ar[i];
			}
			return true;
		} else {
			return true;
		}
	}
	
	/**
	 * checks whether all elements in {@code ar} are the same
	 * @param ar
	 * @return
	 */
	public static boolean areAllEqual(char[] ar) {
		if (ar.length > 1) {
			char elem = ar[0];
			for (int i = 1; i < ar.length; i++) {
				if (ar[i] != elem) {
					return false;
				}
				elem = ar[i];
			}
			return true;
		} else {
			return true;
		}
	}
	
	/**
	 * checks whether all elements in {@code ar} are the same
	 * @param ar
	 * @return
	 */
	public static boolean areAllEqual(boolean[] ar) {
		if (ar.length > 1) {
			boolean elem = ar[0];
			for (int i = 1; i < ar.length; i++) {
				if (ar[i] != elem) {
					return false;
				}
				elem = ar[i];
			}
			return true;
		} else {
			return true;
		}
	}
	
	/**
	 * returns the index of {@code d} found in {@code x}
	 * <br>if {@code x} does not contain {@code d}, then -1 is returned 
	 * @param x
	 * @param d
	 * @return
	 */
	public static int indexAt(final double[] x, final double d) {
		int i = -1;
		for (int j = 0; j < x.length; j++) {
			if (x[j] == d) {
				i = j;
				break;
			}
		}
		return i;
	}
	
	public static <T> T[] removeEmpty(final T[] ar) {
		List<T> list = new ArrayList<T>();
		for (int a = 0; a < ar.length; a++) {
			if (ar[a] != null) {
				list.add(ar[a]);
			}
		}
		return list.toArray(Arrays.copyOf(ar, list.size()));
	}	

	public static <T> T[] removeWith(T[] ar, T removeContent) {
		List<T> list = new ArrayList<T>();
		for (int a = 0; a < ar.length; a++) {
			if (!ar[a].equals(removeContent)) {
				list.add(ar[a]);
			}
		}
		return list.toArray(Arrays.copyOf(ar, list.size()));
	}	
	
}
