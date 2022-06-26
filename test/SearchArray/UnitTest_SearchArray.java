package SearchArray;

import java.util.Arrays;

import org.junit.Test;

import net.sytes.botg.array.ArUtils;
import net.sytes.botg.array.SearchArray;

public class UnitTest_SearchArray {

	double[] ar = {100, Double.NaN, 129812.9, 120.2174914, Double.NaN};
	double[] ar1 = ArUtils.createRandomDoubleArray(1000);
	double[] ar2 = ArUtils.createRandomDoubleArray(998);
	
	@Test
	public void test01() {
		
		System.out.println(Arrays.toString(SearchArray.isNaN(ar)));
		
	}
	
	@Test
	public void test02() {
		
		double[][] dd = SearchArray.separateDataIntoWindows(ar1, 100, true);
		
		System.out.println(Arrays.toString(dd));
		
	}
	
	@Test
	public void test03() {
		
		double[][] dd = SearchArray.separateDataIntoWindows(ar2, 100, false);
		
		System.out.println(Arrays.toString(dd));
		
	}
}
