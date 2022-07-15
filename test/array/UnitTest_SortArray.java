package array;

import java.util.Arrays;

import org.junit.Test;

import net.sytes.botg.array.SortArray;

public class UnitTest_SortArray {

	@Test
	public void test000() {
		byte[] bs = {0, 1, 0, 1, 0, 1, 0, 1};
		
		int n = 10_000;
		long ns = 1_000_000_000;
		long t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			SortArray.oddToEven(bs);
			//int b = bs.length;
		}
		long t2 = System.nanoTime();
		long t_e = t2 - t1;
		double sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
		
		for (int i = 0; i < n; i++) {
			SortArray.oddToEven2(bs);
			//int b = bs.length;
		}
		
		t2 = System.nanoTime();
		t_e = t2 - t1;
		sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
	}
	
	@Test
	public void testOddEven() {
		
		double[] ar1 = {1, 2, 3, 4, 5}; 
		System.out.println(Arrays.toString(SortArray.even(ar1)));
		
		System.out.println(Arrays.toString(SortArray.odd(ar1)));
		
	}
	
}
