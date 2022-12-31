package math;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Random;

import org.junit.Test;
import org.junit.jupiter.api.DisplayName;

import net.sytes.botg.array.ArUtils;
import net.sytes.botg.array.math.Scalar;
import net.sytes.botg.array.math.Vec2Mat;
import net.sytes.botg.array.math.Vec2Scalar;

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
		
		double[][] wAr = Vec2Mat.fixedWindows(ar, 2, true);
		
		System.out.println(Arrays.deepToString(wAr));
		
	}
	
	@Test
	public void test017() {
		
		double[] ar = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
		
		double[][] wAr = Vec2Mat.fixedWindows(ar, 2, false);
		
		System.out.println(Arrays.deepToString(wAr));
		
	}
	
	@Test
	public void test018() {
		int n = 1_000_000;
		
		double[] ar1 = ArUtils.rand(1_000);
		
		long st = System.nanoTime();
		double s = 0;
		for (int i = 0; i < n; i++) {
			s = Vec2Scalar.sum(ar1);
		}
		long et = System.nanoTime();
		System.out.println("Sum: " + s);
		
		double el = et - st;
		double sp = el / n; 
		
		System.out.println("Sum -> Sampling Period per Summation [ns]: " + sp);
		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			s = Vec2Scalar.sum2(ar1);
		}
		et = System.nanoTime();
		System.out.println("Sum: " + s);
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("Sum2 -> Sampling Period per Summation [ns]: " + sp);
		
	}
	
}
