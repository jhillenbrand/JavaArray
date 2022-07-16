package net.sytes.botg.array;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import net.sytes.botg.array.math.Vec2Scalar;

public class SearchArray {
	
	// Suppress default constructor for noninstantiability
	private SearchArray() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * return the indices of the k highest elements from ar[]
	 * @param ar
	 * @param k
	 * @return
	 */
	public static int[] maxInd(double[] ar, int k) {
		int[] sortInds = SortArray.quicksort2(ar);
		int[] maxInds = new int[k];
		int kk = 0;
		for(int i = sortInds.length - 1; i > 0; i--) {
			kk = kk + 1;
			maxInds[kk - 1] = sortInds[i];
			if(kk == k - 1) {
				return maxInds;
			}
		}
		return null;
	}
	
	/**
	 * returns an array with unique elements based on values in {@code x}
	 * <br>faster than unique2 for smaller arrays
	 * @param x
	 * @return
	 */
	public static double[] unique(final double[] x) {
		ArUtils.checkForNull(x);
		ArUtils.checkForEmpty(x);
		int n = x.length;
		int u = 1;
		double[] tmp = new double[n];
		tmp[0] = x[0];
		for (int i = 1; i < n; i++) {
			boolean found = false;
			for (int j = 0; j < u; j++) {
				if (tmp[j] == x[i]) {
					found = true;
					break;
				}
			}
			if (!found) {
				tmp[u] = x[i];
				++u;
			}
		}
		double[] xu = new double[u];
		System.arraycopy(tmp, 0, xu, 0, u);
		return xu;
	}

	
	/**
	 * returns an array with unique elements based on values in {@code x}
	 * @param x
	 * @return
	 */
	public static double[] unique2(final double[] x) {
		ArUtils.checkForNull(x);
		ArUtils.checkForEmpty(x);
		int n = x.length;
		//HashMap<Double, Double> xu = new LinkedHashMap<Double, Double>();
		HashMap<Double, Double> xu = new HashMap<Double, Double>();
		for (int i = 0; i < n; i++) {
			xu.put(x[i], x[i]);
		}
		return ConvertArray.unwrap(xu.values().toArray(new Double[xu.size()]));
	}
	
	public static boolean isArray(Object o) {
		if (o != null && o.getClass().isArray()) {
			return true;
		} else {
			return false;
		}
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
