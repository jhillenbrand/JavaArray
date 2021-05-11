package net.sytes.botg.array;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

public class SortArray {

	
	public static void reverseByteArray(byte[] array) {
		//long start = System.nanoTime();
		int i = 0;
		byte b;
		int n = array.length;
		for (i = 0; i < n / 2; i++) {
			b = array[i];
			array[i] = array[n - i - 1];
			array[n - i - 1] = b;
		}
		//System.out.println(System.nanoTime() - start + " ns elapsed");
	}
	
	// DEBUGGING NOT WORKING DUE TO LIST IMPLEMENTATION
	public static void reverseByteArray2(byte[] array) {
		//long start = System.nanoTime();
		List<byte[]> bL = Arrays.asList(array);
		Collections.reverse(bL);
		//System.out.println(System.nanoTime() - start + " ns elapsed");
	}

	public static byte[] reverseByteArray3(byte[] array) {
		//long start = System.nanoTime();
		int i, j, n;
		n = array.length;
		j = n;
		byte[] newArray = new byte[j];
		for (i = 0; i < n; i++) {
			newArray[j - 1] = array[i];
			j = j - 1;
		}
		return newArray;
		//System.out.println(System.nanoTime() - start + " ns elapsed");
	}
	
	public static byte[] switchOddToEvenBytes(byte[] array) {
		//long start = System.nanoTime();
		byte[] newArray = new byte[array.length];
		if (array.length < 2) {
			throw new RuntimeException("array must contain at least 2 elements");
		}
		for (int i = 1; i < array.length; i += 2) {
			newArray[i - 1] = array[i];
			newArray[i] = array[i - 1];
		}
		return newArray;
		//System.out.println(System.nanoTime() - start + " ns elapsed");
	}
	
	public static void flip(double[] doubleArray) {
		ArrayUtils.reverse(doubleArray);
		return;
	}
	
	public static void flip(int[] intArray) {
		ArrayUtils.reverse(intArray);
		return;
	}
	
	public static void flip(long[] longArray) {
		ArrayUtils.reverse(longArray);
		return;
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
	
	/**
	 * subject to EPL 2.0 in /lic/LICENSE_VOGELLA_QUICKSORT.txt
	 * @see <a href="https://www.vogella.com/tutorials/JavaAlgorithmsQuicksort/article.html">Link</a>
	 * @param ar
	 */
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
	
	/**
	 * subject to EPL 2.0 in /lic/LICENSE_VOGELLA_QUICKSORT.txt
	 * @see <a href="https://www.vogella.com/tutorials/JavaAlgorithmsQuicksort/article.html">Link</a>
	 * @param ar
	 * @return
	 */
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
