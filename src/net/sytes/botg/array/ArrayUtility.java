package net.sytes.botg.array;
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
public class ArrayUtility {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double[] ar = {1,5,3,6,8,0};
		int inds[] = quicksort2(ar);
		System.out.println(Arrays.toString(inds));
	}
	
	public static void reverseByteArray(byte[] array) {
		long start = System.nanoTime();
		int i = 0;
		byte b;
		int n = array.length;
		for (i = 0; i < n / 2; i++) {
			b = array[i];
			array[i] = array[n - i - 1];
			array[n - i - 1] = b;
		}
		System.out.println(System.nanoTime() - start + " ns elapsed");
	}
	
	public static void reverseByteArray2(byte[] array) {
		long start = System.nanoTime();
		Collections.reverse(Arrays.asList(array));
		System.out.println(System.nanoTime() - start + " ns elapsed");
	}

	public static byte[] reverseByteArray3(byte[] array) {
		long start = System.nanoTime();
		int i, j, n;
		n = array.length;
		j = n;
		byte[] newArray = new byte[j];
		for (i = 0; i < n; i++) {
			newArray[j - 1] = array[i];
			j = j - 1;
		}
		System.out.println(System.nanoTime() - start + " ns elapsed");
		return newArray;
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
	
	public static int[] bytesToInt(byte[] bytes, int numOfBytes) {
		long start = System.nanoTime();
		int[] ints = new int[bytes.length / numOfBytes];
		for (int i = 0; i < ints.length; i++) {
			ints[i] = getInt(bytes, numOfBytes * i);
		}
		System.out.println(System.nanoTime() - start + " ns elapsed");
		return ints;
	}
	
	private static int getInt(byte[] arr, int off) {
		return arr[off] << 8 & 0xFF00 | arr[off + 1] & 0xFF;
	}
	
	/**
	 * method returns true/false whether String str is contained in Array ar
	 * @param ar
	 * @param str
	 * @return true/false
	 */
	public static boolean isStringInArray(String[] ar, String str) {
		List<String> list = Arrays.asList(ar);
		if (list.contains(str)) {
			return true;
		} else {
			return false;
		}
	}
	
	public static double[] addValueToDoubleArrayElements(double[] ar, double d) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] + d;
		}
		return ar;
	}
	
	public static double[] multiplyValueToDoubleArrayElements(double[] ar, double d) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] * d;
		}
		return ar;
	}
	
	public static double[] squareDoubleArrayElements(double[] ar) {
		for (int i = 0; i <= ar.length - 1; i++) {
			ar[i] = ar[i] * ar[i];
		}
		return ar;
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
	
	public static double sum(double ar[]) {
		double sum = 0;
		for(double d : ar) {
			sum = sum + d;
		}
		return sum;
	}
	
	public static double median(double ar[]) {
		Arrays.sort(ar);
		double med = 0;
		if (ar.length % 2 == 0) {
		    med = (ar[ar.length / 2] + ar[ar.length / 2 - 1]) / 2;
		} else {
		    med = ar[ar.length / 2];
		}
		return med;
	}
	
	public static double max(double ar[]) {
		double maxVal = 0;
		for(double d : ar) {
			if(maxVal > d) {
				// do nothing
			} else {
				maxVal = d;
			}
		}
		return maxVal;
	}
	
	public static double min(double ar[]) {
		double minVal = 0;
		for(double d : ar) {
			if(minVal > d) {
				minVal = d;
			} else {
				// do nothing
			}
		}
		return minVal;
	}
	
	public static double sumprod(double ar[]) {
		double sumprod = 0;
		for(double d : ar) {
			sumprod = sumprod * d;
		}
		return sumprod;
	}
	
	public static double sumprod(double ar1[], double ar2[]) {
		double sumprod = 0;
		for(int i = 0; i < ar1.length; i++) {
			sumprod = sumprod + ar1[i] * ar2[i];
		}
		return sumprod;
	}
	
	public static double[] max(double ar[], int k) {
		/**
		 * return the k highest elements from ar[]
		 */
		double[] maxVals = new double[k];
		Arrays.sort(ar);
		int kk = 0;
		for(int i = ar.length - 1; i > 0; i--) {
			kk = kk + 1;
			maxVals[kk - 1] = ar[i];
			if(kk == k) {
				return maxVals;
			}
		}
		return maxVals;
	}
	
	public static int[] maxInd(double ar[], int k) {
		/**
		 * return the indices of the k highest elements from ar[]
		 */
		int[] sortInds = quicksort2(ar);
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
	
	public static double[] merge(final double[] ...arrays ) {
	    int size = 0;
	    for (double[] a: arrays) {
	        size += a.length;
	        double[] res = new double[size];
	        int destPos = 0;
	        for ( int i = 0; i < arrays.length; i++ ) {
	            if ( i > 0 ) destPos += arrays[i-1].length;
	            int length = arrays[i].length;
	            System.arraycopy(arrays[i], 0, res, destPos, length);
	        }
	        return res;
	    }
		return null;
	}
	
	public static double[] bubbleSort(double[] ar) {
		for (int a = ar.length - 1; a > 0; a--) {
			for (int b = 0; b < a; b++) {
				if (ar[b] > ar[b + 1]) {
					double temp = ar[b];
					ar[b] = ar[b + 1];
					ar[b + 1] = temp;
				}
			}
		}
		
		return ar;
	}
	
	public static void quicksort(double[] ar) {
        double[] numbers;
        int number;
		// check for empty or null array
        if (ar == null || ar.length == 0){
            return;
        }
        numbers = ar;
        number = ar.length;
        quicksortRec(0, number - 1, numbers);
    }
	
	public static int[] quicksort2(double[] ar) {
        double[] numbers;
        int[] inds = new int[ar.length];
        int number;
		// check for empty or null array
        if (ar == null || ar.length == 0){
            return null;
        }
        numbers = ar;
        number = ar.length;
        for(int i = 0; i < inds.length; i++) {
        	inds[i] = i;
        }
        quicksortRec2(0, number - 1, numbers, inds);
        return inds;
    }

    private static void quicksortRec2(int low, int high, double[] numbers, int[] inds) {
    	int i = low, j = high;
        // Get the pivot element from the middle of the list
        double pivot = numbers[low + (high-low)/2];

        // Divide into two lists
        while (i <= j) {
            // If the current value from the left list is smaller than the pivot
            // element then get the next element from the left list
            while (numbers[i] < pivot) {
                i++;
            }
            // If the current value from the right list is larger than the pivot
            // element then get the next element from the right list
            while (numbers[j] > pivot) {
                j--;
            }

            // If we have found a value in the left list which is larger than
            // the pivot element and if we have found a value in the right list
            // which is smaller than the pivot element then we exchange the
            // values.
            // As we are done we can increase i and j
            if (i <= j) {
            	exchangeQuickSortElements2(i, j, numbers, inds);
                i++;
                j--;
            }
        }
        // Recursion
        if (low < j) {
            quicksortRec2(low, j, numbers, inds);
        }
        if (i < high) {
        	quicksortRec2(i, high, numbers, inds);
        }
	}

	private static void exchangeQuickSortElements2(int i, int j, double[] numbers, int[] inds) {
		// TODO Auto-generated method stub
		double temp = numbers[i];
        numbers[i] = numbers[j];
        numbers[j] = temp;
        int tempInd = inds[i];
        inds[i] = inds[j];
        inds[j] = tempInd;
	}

	private static void quicksortRec(int low, int high, double[] numbers) {
        int i = low, j = high;
        // Get the pivot element from the middle of the list
        double pivot = numbers[low + (high-low)/2];

        // Divide into two lists
        while (i <= j) {
            // If the current value from the left list is smaller than the pivot
            // element then get the next element from the left list
            while (numbers[i] < pivot) {
                i++;
            }
            // If the current value from the right list is larger than the pivot
            // element then get the next element from the right list
            while (numbers[j] > pivot) {
                j--;
            }

            // If we have found a value in the left list which is larger than
            // the pivot element and if we have found a value in the right list
            // which is smaller than the pivot element then we exchange the
            // values.
            // As we are done we can increase i and j
            if (i <= j) {
            	exchangeQuickSortElements(i, j, numbers);
                i++;
                j--;
            }
        }
        // Recursion
        if (low < j) {
            quicksortRec(low, j, numbers);
        }
        if (i < high) {
        	quicksortRec(i, high, numbers);
        }
    }

    private static void exchangeQuickSortElements(int i, int j, double[] numbers) {
        double temp = numbers[i];
        numbers[i] = numbers[j];
        numbers[j] = temp;
    }	
}
