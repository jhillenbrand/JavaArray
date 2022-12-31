package math;

import java.util.Arrays;
import java.util.Random;

import org.junit.Test;

import net.sytes.botg.array.ArUtils;
import net.sytes.botg.array.math.Vec;
import net.sytes.botg.array.math.Vec2Scalar;

public class UnitTest_Vec {

	private double[] ar1 = {12.32, 231234.0, 123123.023, 123123.09, 123, 1231239, 123123.0213, 12356.089, 978997.0324};
	private double[] ar2 = new Random().doubles(1_000_000).toArray();
	
	
	@Test
	public void test000() {
		
		double[] ar = {1.8, 2.0, 8.1, 9.0, 243.2, 2134.123, 23.0, -1.0, 12.9};
		
		double[] newAr = Vec.windowedSpan(ar, 3, true);
		
		System.out.println(Arrays.toString(newAr));
		
	}
	
	@Test
	public void test010() {
		
		double[] ar = {1.8, 2.0, 8.1, 9.0, 243.2, 2134.123, 23.0, -1.0, 12.9, -5.0, 4.5};
		
		double[] newAr = Vec.windowedSpan(ar, 3, false);
		
		System.out.println(Arrays.toString(newAr));
		
	}
	
	@Test
	public void test020() {
		
		double[] x = {22.2,	22.25,	22.3,	22.25,	22.2,	22.3,	22.2,	22.4,	22.4};
		
		double[] xUp = Vec.upsample2(x, 3);
		
		System.out.println(Arrays.toString(xUp));
	}
	
	@Test
	public void test030() {
		
		double[] x = {22.2,	22.25,	22.3,	22.25,	22.2,	22.3,	22.2,	22.4,	22.4};
		
		double[] y = Vec.zeroRange(x, 2, 5);
		
		System.out.println(Arrays.toString(x));
		System.out.println(Arrays.toString(y));
	}
	
	@Test
	public void test040() {
		double[] ar = ArUtils.linspace(0.0, 10.0, 100);
		
		double[] ar2 = Vec.downsampleBruteForce(ar, 100);
		
		System.out.println(Arrays.toString(ar2));
	}
	
	@Test
	public void test041() {
		double[] ar = ArUtils.linspace(0.0, 10.0, 100);
		
		double[] ar2 = Vec.downsampleBruteForce(ar, 50);
		System.out.println(Arrays.toString(ar2));
	}
	
	@Test
	public void test042() {
		double[] ar = ArUtils.linspace(0.0, 10.0, 100);
		
		double[] ar2 = Vec.downsampleBruteForce(ar, 33);
		System.out.println(Arrays.toString(ar2));
	}
	
	@Test
	public void test043() {
		
		int f = 11;
		
		double[] ar = ArUtils.linspace(0.0, 10.0, 100);
		
		double[] ar2 = Vec.downsampleBruteForce(ar, ar.length / f);
		
		System.out.println(Arrays.toString(ar2));
	}
	
	@Test
	public void test044() {
		
		int m = 5;
		
		double[] ar = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
		
		double[] ar2 = Vec.downsampleMax(ar, m);
		
		System.out.println(Arrays.toString(ar2));
	}
	
	@Test
	public void test050() {
		double[] ar = {0.0, 1.0, 2.0, 1.0, 0.0, -1.0, -1.0, 1.0, 0.0, 1.0};
		
		int[] inds = Vec.zeroCrossings(ar);
		
		System.out.println(Arrays.toString(inds));
		
	}
	
	@Test
	public void test051() {
		double[] ar = {0.0, 1.0, 2.0, 1.0, 0.0, -1.0, -1.0, 1.0, 0.0, 1.0};
		
		int[] inds = Vec.zeroCrossings(ar, 3);
		
		System.out.println(Arrays.toString(inds));
		
	}
	
	@Test
	public void test060() {
		
		double[] ar = {0.0, 0.0, -1.0, -0.5, 0.0, 1.5, 2.0, 2.4, 2.3, 2.3, 2.8, 1.0};
		
		int[] extrema = Vec.findLocalExtrema(ar);
		
		System.out.println(Arrays.toString(extrema));
		
	}
	
	@Test
	public void test070() {
		//System.out.println(MathArray.sum2(ar2));
		Vec2Scalar.sum2(ar2);
	}
	
	
	@Test
	public void test071() {
		//System.out.println(MathArray.sum(ar2));
		Vec2Scalar.sum(ar2);
	}
	
}
