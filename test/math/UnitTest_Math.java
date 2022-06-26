package math;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Random;

import org.junit.Test;
import org.junit.jupiter.api.DisplayName;

import net.sytes.botg.array.ArUtils;
import net.sytes.botg.array.math.Vec2Mat;
import net.sytes.botg.array.math.Vec2Scalar;

public class UnitTest_Math {

	private double[] ar1 = {12.32, 231234.0, 123123.023, 123123.09, 123, 1231239, 123123.0213, 12356.089, 978997.0324};
	private double[] ar2 = new Random().doubles(1_000_000).toArray();
	
	@Test
	public void test01() {
		
		assertEquals(ArUtils.closestExponentForBase2(1024), 10);
		
		assertEquals(ArUtils.closestExponentForBase2(2047), 10);
		
	}
	
	@Test
	public void test02() {
		double[] data = ArUtils.createRandomDoubleArray(10000);
		
	}
	
	@Test
	public void test03() {
		//System.out.println(MathArray.sum2(ar2));
		Vec2Scalar.sum2(ar2);
	}
	
	@Test
	public void test04() {
		//System.out.println(MathArray.sum(ar2));
		Vec2Scalar.sum(ar2);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 1")
	public void test05() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isMonotone(ar, true), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 2")
	public void test06() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isMonotone(ar, false), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 3")
	public void test07() {
		double[] ar = {1.0, 1.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isMonotone(ar, true), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 4")
	public void test08() {
		double[] ar = {1.0, 0.99, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 5")
	public void test09() {
		double[] ar = {1.0, 0.99, 0.9, 0.0, -1.1};
		
		assertEquals(ArUtils.isMonotone(ar, false), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 1")
	public void test10() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isStrictMonotone(ar, true), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 2")
	public void test11() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isStrictMonotone(ar, false), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 3")
	public void test12() {
		double[] ar = {1.0, 1.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isStrictMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 4")
	public void test13() {
		double[] ar = {1.0, 0.99, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isStrictMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 5")
	public void test14() {
		double[] ar = {1.0, 0.99, 0.9, 0.0, -1.1};
		
		assertEquals(ArUtils.isStrictMonotone(ar, false), true);
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
		
		double[] ar1 = ArUtils.createRandomDoubleArray(1_000);
		
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
	
	@Test
	public void test019() {
		
		double d = 1.8;
		
		int i = (int) (d / 1);
		
		System.out.println(i);
		
	}
}
