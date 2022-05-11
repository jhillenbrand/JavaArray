package MathArray;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Random;

import org.junit.Test;
import org.junit.jupiter.api.DisplayName;

import net.sytes.botg.array.math.MathArray;
import net.sytes.botg.array.utils.ArrayUtility;

public class UnitTest_MathArray {

	private double[] ar1 = {12.32, 231234.0, 123123.023, 123123.09, 123, 1231239, 123123.0213, 12356.089, 978997.0324};
	private double[] ar2 = new Random().doubles(1_000_000).toArray();
	
	@Test
	public void test01() {
		
		assertEquals(MathArray.closestExponentForBase2(1024), 10);
		
		assertEquals(MathArray.closestExponentForBase2(2047), 10);
		
	}
	
	@Test
	public void test02() {
		double[] data = ArrayUtility.createRandomArray(10000);
		
	}
	
	@Test
	public void test03() {
		//System.out.println(MathArray.sum2(ar2));
		MathArray.sum2(ar2);
	}
	
	@Test
	public void test04() {
		//System.out.println(MathArray.sum(ar2));
		MathArray.sum(ar2);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 1")
	public void test05() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(MathArray.isMonotone(ar, true), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 2")
	public void test06() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(MathArray.isMonotone(ar, false), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 3")
	public void test07() {
		double[] ar = {1.0, 1.0, 3.0, 4.0, 4.1};
		
		assertEquals(MathArray.isMonotone(ar, true), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 4")
	public void test08() {
		double[] ar = {1.0, 0.99, 3.0, 4.0, 4.1};
		
		assertEquals(MathArray.isMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 5")
	public void test09() {
		double[] ar = {1.0, 0.99, 0.9, 0.0, -1.1};
		
		assertEquals(MathArray.isMonotone(ar, false), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 1")
	public void test10() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(MathArray.isStrictMonotone(ar, true), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 2")
	public void test11() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(MathArray.isStrictMonotone(ar, false), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 3")
	public void test12() {
		double[] ar = {1.0, 1.0, 3.0, 4.0, 4.1};
		
		assertEquals(MathArray.isStrictMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 4")
	public void test13() {
		double[] ar = {1.0, 0.99, 3.0, 4.0, 4.1};
		
		assertEquals(MathArray.isStrictMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 5")
	public void test14() {
		double[] ar = {1.0, 0.99, 0.9, 0.0, -1.1};
		
		assertEquals(MathArray.isStrictMonotone(ar, false), true);
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
		
		double[][] wAr = MathArray.fixedWindows(ar, 2, true);
		
		System.out.println(Arrays.deepToString(wAr));
		
	}
	
	@Test
	public void test017() {
		
		double[] ar = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
		
		double[][] wAr = MathArray.fixedWindows(ar, 2, false);
		
		System.out.println(Arrays.deepToString(wAr));
		
	}
}
