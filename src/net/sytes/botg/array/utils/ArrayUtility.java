package net.sytes.botg.array.utils;

import java.util.Arrays;
import java.util.Random;

public class ArrayUtility {

	private static final long BETTER_OF_AS_STREAM_SIZE = 100_000_000;
	
	public static void print2DArray(Object[][] ar) {
		for (int i = 0; i < ar.length; i++) {
			System.out.println(Arrays.toString(ar[i]));
		}
	}
	
	public static void print2DArray(double[][] ar) {
		for (int i = 0; i < ar.length; i++) {
			System.out.println(Arrays.toString(ar[i]));
		}
	}
	
	public static void print2DArray(int[][] ar) {
		for (int i = 0; i < ar.length; i++) {
			System.out.println(Arrays.toString(ar[i]));
		}
	}
	
	public static double[] createRandomDoubleArray(long size) {
		double[] data;
		if (size < BETTER_OF_AS_STREAM_SIZE && size < Integer.MAX_VALUE) {
			data = new double[(int) size];
			Random r = new Random();
			for (int i = 0; i < size; i++) {
				data[i] = r.nextDouble();
			}
		} else {
			data = new Random().doubles(size).toArray();
		}
		return data;
	}
	
	public static int[] createRandomIntArray(long size) {
		int[] data;
		if (size < BETTER_OF_AS_STREAM_SIZE && size < Integer.MAX_VALUE) {
			data = new int[(int) size];
			Random r = new Random();
			for (int i = 0; i < size; i++) {
				data[i] = r.nextInt();
			}
		} else {
			data = new Random().ints(size).toArray();
		}
		return data;
	}
	
	/**
	 * returns a double array starting from {@code start} to {@code end} with {@code step}
	 * @param start
	 * @param end
	 * @param step
	 * @return
	 */
	public static double[] linspace(double start, double end, double step) {
		int size = (int) ((end - start ) / step);
		double[] data = new double[size];
		for (int i = 0; i < size; i++) {
			data[i] = start + step * i;
		}
		return data;
	}
	
	/**
	 * returns a double array starting from {@code start} to {@code end} with equally spaced steps, so that the {@code size} is met
	 * @param start
	 * @param end
	 * @param size
	 * @return
	 */
	public static double[] linspace(double start, double end, int size) {
		double[] data = new double[size];
		double step = (end - start) / (size - 1);
		for (int i = 0; i < size; i++) {
			data[i] = start + step * i;
		}
		return data;
	}
	
	/**
	 * returns a double array starting from 0 and incrementing by 1 for {@code size} times
	 * @param size
	 * @return
	 */
	public static double[] linspace(int size) {
		return linspace(0, size, size);
	}
	
	public static double[] nan(int size) {
		double[] data = new double[size];
		for (int i = 0; i < size; i++) {
			data[i] = Double.NaN;
		}
		return data;
	}	

	/**
	 * add the element x to the end of the array ar and return the new array
	 * @param <T>
	 * @param x
	 * @param ar
	 * @return
	 */
	public static <T> T[] add(T x, T[] ar) {
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
	public static double[] add(double d, double[] ar) {
		double[] newAr = new double[ar.length + 1];
		System.arraycopy(ar, 0, newAr, 0, ar.length);
		newAr[ar.length] = d;
		return newAr;
	}
}
