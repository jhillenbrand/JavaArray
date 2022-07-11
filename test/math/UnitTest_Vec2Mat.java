package math;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.ArUtils;
import net.sytes.botg.array.math.Vec2Mat;

public class UnitTest_Vec2Mat {

	@Test
	public void test000() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
		int numSamples = x.length;
		int fftSampleSize = 4;
		int overlapFactor = 2;
		int a = 0;
		
		if (overlapFactor > 1) {

			int numOverlappedSamples = numSamples * overlapFactor;
			int backSamples = fftSampleSize * (overlapFactor - 1) / overlapFactor;
			int fftSampleSize_1 = fftSampleSize - 1;
			double[] overlapAmp = new double[numOverlappedSamples];
			a = 0;
			for (int i = 0; i < numSamples; i++) {
				overlapAmp[a++] = x[i];
				if (a % fftSampleSize == fftSampleSize_1) {
					// overlap
					i -= backSamples;
				}
			}
			numSamples = numOverlappedSamples;
			x = overlapAmp;
		}
		
		System.out.println(Arrays.toString(x));
		
	}
	
	@Test
	public void test010() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
		
		double[][] w = Vec2Mat.overlapWindows(x, 4, 2);
		
		ArUtils.print(w);
		
	}	
	
	@Test
	public void test011() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] w = Vec2Mat.overlapWindows(x, 4, 2);
		
		ArUtils.print(w);
		
	}	
	
	@Test
	public void test012() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] w = Vec2Mat.overlapWindows(x, 4, 1);
		
		ArUtils.print(w);
		
	}	
	
	@Test
	public void test013() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] w = Vec2Mat.overlapWindows(x, 3, 2);
		
		ArUtils.print(w);
		
	}	
	
	@Test
	public void test014() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] w = Vec2Mat.overlapWindows(x, 3, 1);
		
		ArUtils.print(w);
		
	}	
	
	@Test
	public void test015() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] w = Vec2Mat.overlapWindows(x, 5, 2);
		
		ArUtils.print(w);
		
	}	
	
	@Test
	public void test016() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] w = Vec2Mat.overlapWindows(x, 5, 3);
		
		ArUtils.print(w);
		
	}	
	
	@Test
	public void test017() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] w = Vec2Mat.overlapWindows(x, 5, 0);
		
		ArUtils.print(w);
		
	}	
	
	@Test
	public void test018() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] w = Vec2Mat.overlapWindows(x, 5, 4);
		
		ArUtils.print(w);
		
	}	
	
	@Test
	public void test019() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
		
		double[][] w = Vec2Mat.overlapWindows(x, 3, 2);
		
		ArUtils.print(w);
		
	}	
}
