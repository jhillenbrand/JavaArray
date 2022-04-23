package MathArray;

import java.util.Random;

import org.junit.Test;

import net.sytes.botg.array.math.MathArray;

public class UnitTest_ScalarProd {
	
	private double[] ar1 = new Random().doubles(10_000_000).toArray();
	private double[] ar2 = new Random().doubles(10_000_000).toArray();
	
	@Test
	public void testScalarProd00() {
		System.out.println("WARMUP");
	}
	
	@Test
	public void testScalarProd01() {
		long t1 = System.nanoTime();
		MathArray.scalarProd(ar1, ar2);
		long t2 = System.nanoTime();
		System.out.println("01 - " + (t2 - t1) + " ns");
	}
}
