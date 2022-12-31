package array;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.Test;
import org.junit.jupiter.api.DisplayName;

import net.sytes.botg.array.ArUtils;
import net.sytes.botg.array.ConvertArray;
import net.sytes.botg.array.math.Vec;

public class UnitTest_ArUtils {

	private static final int SIZE = 10_000_000;
	
	@Test
	public void test00() {
		
		System.out.println("WARMUP");
		System.out.println(Integer.MAX_VALUE);
	}	
	

	@Test
	public void test010() {
		double[] data = ArUtils.rand(10000);
		
	}
	
	@Test
	public void test011() {
		
		ArUtils.rand(SIZE);
		
	}
	
	@Test
	public void test020() {
		
		Double[] ar = {10.0, 128.234, 3984.123};
		Double d = 10.9877;
		
		System.out.println(Arrays.toString(ar));
		System.out.println(d);
		
		Double[] nAr = Vec.append(ar, d);
		System.out.println(Arrays.toString(nAr));
	}
	
	@Test
	public void test030() {
		
		double[] ar = {10.0, 128.234, 3984.123};
		double d = 10.9877;
		
		System.out.println(Arrays.toString(ar));
		System.out.println(d);
		
		double[] nAr = Vec.append(ar, d);
		System.out.println(Arrays.toString(nAr));
	}
	
	@Test
	public void test040() {
		
		Character[] cw = new Character[2];
		cw[0] = 'A';
		cw[1] = 'B';
		
		char[] c = ConvertArray.unwrap(cw);
		
		System.out.println(c);		
	}
	
	@Test
	public void test050() {
		
		double[] d = ArUtils.rand(100);
		
		double[] dd = ArUtils.subArray(d, 0, 100);
		
		System.out.println(Arrays.toString(dd));
		
	}
	
	@Test
	public void test060() {
		
		double[] d = ArUtils.rand(100);
		
		double[] dd = ArUtils.subArray(d, 1, 90);
		
		System.out.println(Arrays.toString(dd));
		
	}
	
	@Test
	public void test070() {
		
		double[] x = ArUtils.linspace(-1.0, 1.0, 10);
		double[] y = ArUtils.linspace(-1.0, 1.0, 5);
		
		double[][][] XY = ArUtils.meshgrid(x, y);
				
		ArUtils.print(XY[0]);
		ArUtils.print(XY[1]);
		
	}
	
	@Test
	public void test080() {
		
		int s = 1_000;
		int n = 100_000;		

		//int[] ints = ArUtils.createRandomIntArray(s);
		int[] ints = ArUtils.linspace(1, 100);
				
		long t1 = System.nanoTime();
		
		
		for (int i = 0; i < ints.length; i++) {
			ArUtils.checkForGreaterEqualZero2(ints);
		}
		
		long t2 = System.nanoTime();
		long t_e = t2 - t1;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		
		
		
		t1 = System.nanoTime();		
		for (int i = 0; i < ints.length; i++) {
			ArUtils.checkForGreaterEqualZero(ints);
		}
		
		t2 = System.nanoTime();
		t_e = t2 - t1;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
	}

	@Test
	public void test090() {
		double[] nans = ArUtils.nan(10);
		
		System.out.println(Arrays.toString(nans));
		
	}
	
	@Test
	public void test091() {
		double[][] nans = ArUtils.nan(10, 2);
		
		System.out.println(Arrays.deepToString(nans));
		
	}
	
	@Test
	public void test100() {
		double[] nans = ArUtils.nan(10);
		double[] rands = ArUtils.rand(10);
		
		double[] prods = Vec.prod(rands, nans);
		
		System.out.println(Arrays.toString(prods));
		
	}
	
	@Test
	public void test110() {
		
		double[] logspace = ArUtils.logspace(1, 10, 10);
		System.out.println(Arrays.toString(logspace));
		
	}
	

	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 1")
	public void test120() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isMonotone(ar, true), true);
	}
	

	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 2")
	public void test121() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isMonotone(ar, false), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 3")
	public void test122() {
		double[] ar = {1.0, 1.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isMonotone(ar, true), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 4")
	public void test123() {
		double[] ar = {1.0, 0.99, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 5")
	public void test124() {
		double[] ar = {1.0, 0.99, 0.9, 0.0, -1.1};
		
		assertEquals(ArUtils.isMonotone(ar, false), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 1")
	public void test125() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isStrictMonotone(ar, true), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 2")
	public void test126() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isStrictMonotone(ar, false), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 3")
	public void test127() {
		double[] ar = {1.0, 1.0, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isStrictMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 4")
	public void test128() {
		double[] ar = {1.0, 0.99, 3.0, 4.0, 4.1};
		
		assertEquals(ArUtils.isStrictMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 5")
	public void test129() {
		double[] ar = {1.0, 0.99, 0.9, 0.0, -1.1};
		
		assertEquals(ArUtils.isStrictMonotone(ar, false), true);
	}
	
	@Test
	public void test130() {
		
		int n = 100;
		int r = 20_000;
		
		double[][] X = ArUtils.unitMatrix(n);
		double[][] Y = null;
		
		long st = System.nanoTime();
		for (int i = 0; i < r; i++) {
			
			Y = ArUtils.copy2(X);
			
		}
		
		long et = System.nanoTime();
		double el = et - st;
		double sp = el / n; 
		
		System.out.println("copy2 -> Sampling Period per Element [ns]: " + sp);
		
		st = System.nanoTime();
		for (int i = 0; i < r; i++) {
			
			Y = ArUtils.copy(X);
			
		}
		
		et = System.nanoTime();
		el = et - st;
		sp = el / n; 
		
		System.out.println("copy -> Sampling Period per Element [ns]: " + sp);
		
		System.out.println("X = [" + X[0].length + ", " + X.length + "], Y = [" + Y[0].length + ", " + Y.length + "]");
		
	}
	
	@Test
	public void test140() {
		
		double[] x = new double[] {1.0, 2.0, 3.0};
		
		double[] y = ArUtils.copy(x);
		
		y[0] = 0.0;
		
		System.out.println(Arrays.toString(x));
		
	}
	
	@Test
	public void test150() {
		double[][] X = ArUtils.incrementMat(3,  3);
		ArUtils.print(X);
	}
	
	@Test
	public void test151() {
		double[][] X = ArUtils.incrementMat(2,  5);
		ArUtils.print(X);
	}
	
	@Test
	public void test152() {
		double[][] X = ArUtils.incrementMat(5,  3);
		ArUtils.print(X);
	}
	
}
