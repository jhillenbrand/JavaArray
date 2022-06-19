package net.sytes.botg.array;

import java.util.Arrays;
import java.util.List;

import net.sytes.botg.array.math.Vec2Scalar;

public class SearchArray {
	
	// Suppress default constructor for noninstantiability
	private SearchArray() {
		throw new AssertionError();
	}
	
	public static boolean isArray(Object o) {
		if (o != null && o.getClass().isArray()) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * returns the array {@code data} split into windows of size {@code window}
	 * if dropUneven is set to {@code false}, then the remaining elements are filled with NaN
	 * @param data
	 * @param window
	 * @param dropUneven
	 * @return double[][]
	 */
	public static double[][] separateDataIntoWindows(double[] data, int window, boolean dropUneven){
		double[][] dataWindows = null;
		if (dropUneven) {
			int w = data.length / window;
			dataWindows = new double[w][window];
			int s = 0;
			int e = window;
			for (int i = 0; i < w; i++) {
				s = window * i;
				System.arraycopy(data, s, dataWindows[i], 0, e);
			}
		} else {
			int w = data.length / window;
			int rem = data.length % window;
			if (rem > 0) {
				 ++w;
			}
			dataWindows = new double[w][window];
			int s = 0;
			int e = window -1;
			for (int i = 0; i < w; i++) {
				if (i == w - 1 && rem > 0) {
					double[] nanAr = ArrayUtility.nan(rem);
					s = window * i;
					System.arraycopy(data, s, dataWindows[i], 0, rem);
					System.arraycopy(nanAr, 0, dataWindows[i], rem, e - rem + 1);
				} else {
					s = window * i;
					System.arraycopy(data, s, dataWindows[i], 0, e);
				}
			}			
		}
		return dataWindows;
	}
	
	/**
	 * returns a boolean[] array containing true for each element in {@code data} that is within the bounds [{@code lowerLimit}, {@code upperLimit}]
	 * @param data
	 * @param lowerLimit
	 * @param upperLimit
	 * @return
	 */
	public static boolean[] isInRange(double[] data, double lowerLimit, double upperLimit) {
		boolean[] inds = new boolean[data.length];
		for (int i = 0; i < data.length; i++) {
			if (data[i] >= lowerLimit && data[i] <= upperLimit) {
				inds[i] = true;
			}
		}
		return inds;
	}
	
	/**
	 * returns the elements specified with {@code inds} in {@code data} into new array
	 * @param inds
	 * @param data
	 * @return double[]
	 */
	public static double[] elementsAt(boolean[] inds, double[] data) {
		int n = Vec2Scalar.sum(inds);
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
	 * returns an boolean[] array that is true for every element in {@code data} with {@code NaN}
	 * @param data
	 * @return boolean[]
	 */
	public static boolean[] isNaN(double[] data) {
		boolean[] inds = new boolean[data.length];
		for (int i = 0; i < data.length; i++) {
			// returns true only if the element is NaN
			if (data[i] != data[i]) {
				inds[i] = true;
			}
		}
		return inds;
	}
	
	/**
	 * method returns true/false whether String str is contained in Array ar
	 * @param ar
	 * @param str
	 * @return true/false
	 */
	public static boolean isStringInArray(String[] ar, String str) {
		// TODO this is slow, because contains also iterates over list
		List<String> list = Arrays.asList(ar);
		if (list.contains(str)) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * method returns true/false whether any of the array elements ar are contained in String s
	 * @param s
	 * @param ar
	 * @return 
	 */
	public static boolean containsStringArrayElement(String s, String[] ar) {
		for (String a : ar) {
			if (s.contains(a)) {
				return true;
			}
		}
		return false;
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
