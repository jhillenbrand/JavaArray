package testing.scalarprod;

import java.util.Random;

import org.junit.Test;

import net.sytes.botg.array.MathArray;

public class UnitTest_ScalarProd {
	
	private double[] ar1 = new Random().doubles(10_000_000).toArray();
	private double[] ar2 = new Random().doubles(10_000_000).toArray();
	
	@Test
	public void testScalarProd() {
		MathArray.scalarProd(ar1, ar2);
	}
	
	@Test
	public void testScalarProd2() {
		MathArray.scalarProd2(ar1, ar2);
	}
}
