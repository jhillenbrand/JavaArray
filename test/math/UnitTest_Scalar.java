package math;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.Test;

import net.sytes.botg.array.math.Scalar;
import net.sytes.botg.array.math.Vec;

public class UnitTest_Scalar {

	@Test
	public void test01() {
		
		assertEquals(Scalar.closestExponentForBase2(1024), 10);
		
		assertEquals(Scalar.closestExponentForBase2(2047), 10);
		
	}
	
	@Test
	public void test015() {
		
		int i = 3;
		int j = 2;
		
		System.out.println("" + (i / j));
	}
	
	@Test
	public void test016() {
		
		double[] ar = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
		
		double[][] wAr = Vec.fixedWindows(ar, 2, true);
		
		System.out.println(Arrays.deepToString(wAr));
		
	}
	
	@Test
	public void test017() {
		
		double[] ar = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
		
		double[][] wAr = Vec.fixedWindows(ar, 2, false);
		
		System.out.println(Arrays.deepToString(wAr));
		
	}
	
	@Test
	public void test018() {
		int n = 1_000_000;
		
		double[] ar1 = Vec.rand(1_000);
		
		long st = System.nanoTime();
		double s = 0;
		for (int i = 0; i < n; i++) {
			s = Vec.sum(ar1);
		}
		long et = System.nanoTime();
		System.out.println("Sum: " + s);
		
		double el = et - st;
		double sp = el / n; 
		
		System.out.println("Sum -> Sampling Period per Summation [ns]: " + sp);
		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			s = Vec.sum2(ar1);
		}
		et = System.nanoTime();
		System.out.println("Sum: " + s);
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("Sum2 -> Sampling Period per Summation [ns]: " + sp);
		
	}
	
	@Test
	public void test020() {
		
		double x1 = 1.4;
		double x2 = 1.5;
		
		System.out.println("round(" + x1 + ") --> " + Scalar.round(x1));
		System.out.println("round(" + x2 + ") --> " + Scalar.round(x2));
		
	}
	
}
