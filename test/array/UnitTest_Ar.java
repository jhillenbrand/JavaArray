package array;

import org.junit.Test;
import org.junit.jupiter.api.Assertions;

import net.sytes.botg.array.Ar;
import net.sytes.botg.array.math.Vec;

public class UnitTest_Ar {


	
	@Test
	public void test000() {
		int[] d = Vec.randInt(1_000_000);		
		Ar.intToDouble2(d);
	}
	
	@Test
	public void test001() {
		int[] d = Vec.randInt(1_000_000);		
		Ar.toDouble(d);
	}
	
	@Test
	public void test002() {
		int[] d = Vec.randInt(10_000_000);		
		Ar.intToDouble2(d);
	}
	
	@Test
	public void test003() {
		int[] d = Vec.randInt(10_000_000);		
		Ar.toDouble(d);
	}
	
	@Test
	public void test004() {
		int[] ints = Vec.randInt(100_000);
		
		int n = 10_000;
		long ns = 1_000_000_000;
		long t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			double[] bs = Ar.toDouble(ints);
			//int b = bs.length;
		}
		long t2 = System.nanoTime();
		long t_e = t2 - t1;
		double sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
		
		for (int i = 0; i < n; i++) {
			double[] bs = Ar.intToDouble2(ints);
			//int b = bs.length;
		}
		
		t2 = System.nanoTime();
		t_e = t2 - t1;
		sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
	}
	
	@Test
	public void test010() {
		final int len = 100;
		Object[] objAr = new Object[len]; 
		for (int i = 0; i < len; i++) {
			objAr[i] = i * 0.5;
		}
		
		int n = 10_000;
		long ns = 1_000_000_000;
		long t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			double[] dAr = Ar.toDouble(objAr);
			//int b = bs.length;
		}
		long t2 = System.nanoTime();
		long t_e = t2 - t1;
		double sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
		
		t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			double[] dAr = Ar.toDouble2(objAr);
			//int b = bs.length;
		}
		t2 = System.nanoTime();
		t_e = t2 - t1;
		sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
	}
	
	@Test
	public void test020() {
		
		Character[] cw = new Character[2];
		cw[0] = 'A';
		cw[1] = 'B';
		
		char[] c = Ar.unwrap(cw);
		
		System.out.println(c);		
	}
	
	@Test
	public void test030() {
		
		double[] x = Vec.linspace(10);
		
		double[] x_ = Ar.removeAt(x, 1);
		
		Ar.print(x);
		Ar.print(x_);
	}
	
	@Test
	public void test050() {
		
		double[] x = new double[] {0.0, 1.0, 9.0, 1.0, 5.0, 4.0, 1.0};
		
		double d = 1.0;
		
		int[] inds = Ar.find(x, d);
		
		Ar.print(inds);
		
	}
	
	@Test
	public void test060() {
		
		double[][] X = {{1, 1, 1},
						{1, 0, 1},
						{0, 1, 1},
						{1, 1, 1}};
		
		double[][] Y = Ar.unique(X);
		
		Ar.print(Y);
		
	}
	
	@Test
	public void test070() {
		
		String[] ar1 = {"a", "A", "a"}; 
		
		Assertions.assertFalse(Ar.areAllEqual(ar1));
		
		String[] ar2 = {"a", "a", "a"}; 
		
		Assertions.assertTrue(Ar.areAllEqual(ar2));
		
		Integer[] ar3 = {1, 1, 1};
		
		Assertions.assertTrue(Ar.areAllEqual(ar3));
		
		int[] ar4 = {1, 1, 1};
		
		Assertions.assertTrue(Ar.areAllEqual(ar4));
		
		
	}
	
}
