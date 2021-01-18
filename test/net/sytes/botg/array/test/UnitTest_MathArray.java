package net.sytes.botg.array.test;

import java.util.Random;

import org.junit.Test;

import net.sytes.botg.array.MathArray;

public class UnitTest_MathArray {

	private double[] ar1 = {12.32, 231234.0, 123123.023, 123123.09, 123, 1231239, 123123.0213, 12356.089, 978997.0324};
	
	private double[] ar2 = new Random().doubles(10_000_000).toArray();

	@Test
	public void testSummation1() {
		System.out.println(MathArray.sum2(ar2));
	}
	
	@Test
	public void testSummation2() {
		System.out.println(MathArray.sum(ar2));
	}
	
	
}
