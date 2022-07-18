package array;

import java.util.Arrays;

import org.junit.Test;

import net.sytes.botg.array.ArUtils;
import net.sytes.botg.array.SearchArray;
import net.sytes.botg.array.math.Vec2Mat;

public class UnitTest_SearchArray {

	double[] ar = {100, Double.NaN, 129812.9, 120.2174914, Double.NaN};
	double[] ar1 = ArUtils.rand(1000);
	double[] ar2 = ArUtils.rand(998);
	
	@Test
	public void test01() {
		
		System.out.println(Arrays.toString(SearchArray.isNaN(ar)));
		
	}
	
	@Test
	public void test02() {
		
		double[][] dd = Vec2Mat.separateDataIntoWindows(ar1, 100, true);
		
		System.out.println(Arrays.toString(dd));
		
	}
	
	@Test
	public void test03() {
		
		double[][] dd = Vec2Mat.separateDataIntoWindows(ar2, 100, false);
		
		System.out.println(Arrays.toString(dd));
		
	}
	
	@Test
	public void test04() {
		
		double[] x = {1.0, 1.1, 9.1, 76.3, 1.1, 2.0, 9.1};
		
		double[] xu = SearchArray.unique(x);
		
		System.out.println(Arrays.toString(xu));
		
	}
	
	@Test
	public void test05() {
		
		double[] x = {1.0, 1.1, 9.1, 76.3, 1.1, 2.0, 9.1};
		
		double[] xu = SearchArray.unique2(x);
		
		System.out.println(Arrays.toString(xu));
		
	}
	
	@Test
	public void test06() {
		
		double[] x = {1.0, 1.1, 9.1, 76.3, 1.1, 2.0, 9.1};
		
		int n = 1_000_000;
		long t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			double[] xu = SearchArray.unique2(x);
			//int b = bs.length;
		}
		long t2 = System.nanoTime();
		long t_e = t2 - t1;
		
		System.out.println("unique2: " + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		
		
		t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			double[] xu = SearchArray.unique(x);
			//int b = bs.length;
		}
		t2 = System.nanoTime();
		t_e = t2 - t1;
		
		System.out.println("unique: " + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
				
	}
}
