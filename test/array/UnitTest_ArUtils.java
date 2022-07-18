package array;

import java.util.Arrays;

import org.junit.Test;

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
	public void test01() {
		
		ArUtils.rand(SIZE);
		
	}
	
	@Test
	public void test02() {
		
		Double[] ar = {10.0, 128.234, 3984.123};
		Double d = 10.9877;
		
		System.out.println(Arrays.toString(ar));
		System.out.println(d);
		
		Double[] nAr = Vec.append(ar, d);
		System.out.println(Arrays.toString(nAr));
	}
	
	@Test
	public void test03() {
		
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
}
