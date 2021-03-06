package net.sytes.botg.array;
import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * 
 */

/**
 * @author jonas
 *
 */
public class ConvertArray {
	
	public static byte[] convertDoubleArrayToByteArray(double[] ar, int doublePrecision) {
		byte[] byteArray = new byte[ar.length * doublePrecision / 8];
		for (int a = 0; a < ar.length; a++) {
			byte[] bytes = convertDoubleToBytes(ar[a], doublePrecision);
			for (int b = 0; b < bytes.length; b++) {
				byteArray[a * doublePrecision + b] = bytes[b];
			}
		}
		return byteArray;
	}
	
	public static byte[] convertNumberArrayToByteArray(Number[] ar, int doublePrecision) {
		byte[] byteArray = new byte[ar.length * doublePrecision / 8];
		int index = 0;
		int numOfNumbers = ar.length;
		for (int a = 0; a < numOfNumbers; a++) {
			byte[] bytes = convertNumberToBytes(ar[a], doublePrecision);
			int numOfBytes = bytes.length;
			//if (a == numOfNumbers - 1) {
			//	System.out.println("DEBUG NOTE");
			//}
			for (int b = 0; b < numOfBytes; b++) {
				index = a * doublePrecision / 8 + b;
				byteArray[index] = bytes[b];
			}
		}
		return byteArray;
	}
	
	/**
	 * method returns the Number in byte[], assuming 64bit double precision
	 * @param d
	 * @return
	 */
	public static byte[] convertNumberToBytes(Number d) {
		return convertNumberToBytes(d.doubleValue(), 64);
	}
	
	/**
	 * method returns the Number in byte[], assuming <doublePrecision [int]>-bit precision
	 * @param d
	 * @param doublePrecision
	 * @return
	 */
	public static byte[] convertNumberToBytes(Number d, int doublePrecision) {
		return convertDoubleToBytes2(d.doubleValue(), doublePrecision);
	}
	
	/**
	 * 
	 * @param d
	 * @param doublePrecision in bits (e.g. 64); 
	 * @return
	 */
	public static byte[] convertDoubleToBytes(double d, int doublePrecision) {
		byte[] bytes = new byte[doublePrecision / 8];
		long lng = Double.doubleToLongBits(d);
		for(int i = 0; i < bytes.length; i++) {
			bytes[i] = (byte)((lng >> ((7 - i) * 8)) & 0xff);
		}
		return bytes;
	}
	
	public static byte[] convertDoubleToBytes(double d, int size, double precision) {
		double temp = d;
		temp = temp * Math.pow(10, precision);
		/*
		for (int i = 0; i < precision; i++) {
			temp *= 10;
		}
		*/
	    int output = (int) temp;
	    String strOut = String.format("%0" + size + "d", output);
	    return strOut.getBytes();
	}

	/**
	 * 
	 * @param d
	 * @param precision
	 * @return
	 */
	public static byte[] convertDoubleToBytes2(double d, int precision) {
		byte[] bytes = new byte[precision / 8];
	    ByteBuffer.wrap(bytes).putDouble(d);
	    return bytes;
	}
	
	/**
	 * 
	 * @param bytes
	 * @return
	 */
	public static double convertBytesToDouble(byte[] bytes) {
	    return ByteBuffer.wrap(bytes).getDouble();
	}
	
	/**
	 * 
	 * @param s
	 * @return
	 */
	public static byte[] convertShortToBytes(short s) {
		byte[] bytes = new byte[2];
		ByteBuffer.wrap(bytes).putShort(s);
		return bytes;
	}
		
	/**
	 * 
	 * @param i
	 * @param precision
	 * @return
	 */
	public static byte[] convertIntegerToBytes(int i, int precision) {
		BigInteger bigInt = BigInteger.valueOf(i);
		return bigInt.toByteArray();
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
	
	public static char[] bytesToHex(byte[] bytes) {
		long start = System.nanoTime();
	    char[] hexArray = "0123456789ABCDEF".toCharArray();
	    char[] hexChars = new char[bytes.length * 2];
	    for ( int j = 0; j < bytes.length; j++ ) {
	        int v = bytes[j] & 0xFF;
	        hexChars[j * 2] = hexArray[v >>> 4];
	        hexChars[j * 2 + 1] = hexArray[v & 0x0F];
	    }
	    System.out.println(System.nanoTime() - start + " ns elapsed");
	    return hexChars;
	}
	
	public static int[] bytesToInt(byte[] bytes, int precision) {
		//long start = System.nanoTime();
		int[] ints = new int[bytes.length / precision];
		for (int i = 0; i < ints.length; i++) {
			ints[i] = getInt(bytes, precision * i);
		}
		//System.out.println(System.nanoTime() - start + " ns elapsed");
		return ints;
	}
	
	public static double bytesArrayToDouble(byte[] bytes) {
		return ByteBuffer.wrap(bytes).getDouble();
	}
	
	/**
	 * assuming 64bit double precision a double array is returned based on element wise conversion of 8 byte of byteArray
	 * @param byteArray
	 * @return double[]
	 */
	public static double[] byteArrayToDoubleArray(byte[] byteArray) {
		double[] doubleArray = new double[byteArray.length / 8]; // divided by 8, because 8 * 8 = 64 bit precision
		byte[] bytes = null;		
		for (int i = 0; i < doubleArray.length; i++) {
			bytes = Arrays.copyOfRange(byteArray, i * 8, (i + 1) * 8);
			doubleArray[i] = bytesArrayToDouble(bytes);
		}
		return doubleArray;		
	}
	
	private static int getInt(byte[] arr, int off) {
		return arr[off] << 8 & 0xFF00 | arr[off + 1] & 0xFF;
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
	
    public static int countStringOccurencesInArray(String[] ar, String s, boolean caseSensitive) {
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
    
    public static boolean isIntInArray(int[] ar, int i) {
    	for (int j = 0; j < ar.length; j++) {
    		if (ar[j] == i) {
    			return true;
    		}
    	}
    	return false;
    }
}
