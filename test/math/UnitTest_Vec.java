package math;

import java.util.Arrays;

import org.junit.Test;

import net.sytes.botg.array.math.Vec;
import net.sytes.botg.array.math.Vec2Mat;

public class UnitTest_Vec {

	
	@Test
	public void test000() {
		
		double[] ar = {1.8, 2.0, 8.1, 9.0, 243.2, 2134.123, 23.0, -1.0, 12.9};
		
		double[] newAr = Vec.windowedSpan(ar, 3, true);
		
		System.out.println(Arrays.toString(newAr));
		
	}
	
	@Test
	public void test010() {
		
		double[] ar = {1.8, 2.0, 8.1, 9.0, 243.2, 2134.123, 23.0, -1.0, 12.9, -5.0, 4.5};
		
		double[] newAr = Vec.windowedSpan(ar, 3, false);
		
		System.out.println(Arrays.toString(newAr));
		
	}
	
	@Test
	public void test020() {
		
		double[] x = {22.2,	22.25,	22.3,	22.25,	22.2,	22.3,	22.2,	22.4,	22.4};
		
		double[] xUp = Vec.upsample2(x, 3);
		
		System.out.println(Arrays.toString(xUp));
	}
	
	@Test
	public void test030() {
		
		double[] x = {22.2,	22.25,	22.3,	22.25,	22.2,	22.3,	22.2,	22.4,	22.4};
		
		double[] y = Vec.zeroRange(x, 2, 5);
		
		System.out.println(Arrays.toString(x));
		System.out.println(Arrays.toString(y));
	}
	
}