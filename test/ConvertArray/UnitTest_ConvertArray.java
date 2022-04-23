package ConvertArray;

import org.junit.Test;

import net.sytes.botg.array.ConvertArray;
import net.sytes.botg.array.utils.ArrayUtility;

public class UnitTest_ConvertArray {

	@Test
	public void test000() {
		int[] d = ArrayUtility.createRandomIntArray(1_000_000);		
		ConvertArray.intToDouble2(d);
	}
	
	@Test
	public void test001() {
		int[] d = ArrayUtility.createRandomIntArray(1_000_000);		
		ConvertArray.intToDouble(d);
	}
	
	@Test
	public void test002() {
		int[] d = ArrayUtility.createRandomIntArray(10_000_000);		
		ConvertArray.intToDouble2(d);
	}
	
	@Test
	public void test003() {
		int[] d = ArrayUtility.createRandomIntArray(10_000_000);		
		ConvertArray.intToDouble(d);
	}
	
	@Test
	public void test004() {
		int[] ints = ArrayUtility.createRandomIntArray(100_000);
		
		int n = 10_000;
		long ns = 1_000_000_000;
		long t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			double[] bs = ConvertArray.intToDouble(ints);
			//int b = bs.length;
		}
		long t2 = System.nanoTime();
		long t_e = t2 - t1;
		double sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
		
		for (int i = 0; i < n; i++) {
			double[] bs = ConvertArray.intToDouble2(ints);
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
			double[] dAr = ConvertArray.objectToDouble(objAr);
			//int b = bs.length;
		}
		long t2 = System.nanoTime();
		long t_e = t2 - t1;
		double sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
		
		t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			double[] dAr = ConvertArray.objectToDouble2(objAr);
			//int b = bs.length;
		}
		t2 = System.nanoTime();
		t_e = t2 - t1;
		sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
	}
	
}
