package MathArray;

import java.util.Arrays;

import org.junit.Test;

import net.sytes.botg.array.SortArray;

public class UnitTest_SortArray {

	private double[] ar1 = {1, 2, 3, 4, 5}; 
	
	@Test
	public void testOddEven() {
		
		System.out.println(Arrays.toString(SortArray.even(ar1)));
		
		System.out.println(Arrays.toString(SortArray.odd(ar1)));
		
	}
	
}
