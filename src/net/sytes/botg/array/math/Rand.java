package net.sytes.botg.array.math;

import java.util.Random;

public class Rand {
	
	public static int getRandomNumberInRange(int min, int max) {
		if (min >= max) {
			throw new IllegalArgumentException("max must be greater than min");
		}
		Random r = new Random();
		return r.nextInt((max - min) + 1) + min;
	}
	
	public static int[] getRandomNumberInRange(int n, int min, int max) {
		if (min >= max) {
			throw new IllegalArgumentException("max must be greater than min");
		}
		int[] rInts = new int[n];		
		Random r = new Random();
		for (int i = 0; i < n; i++) {
			rInts[i] = r.nextInt((max - min) + 1) + min;
		}
		return rInts;
	}
	
	public static double getRandomDoubleInRange(double min, double max) {
		Random r = new Random();
		double randomValue = min + (max - min) * r.nextDouble();
		return randomValue;
	}
}
