package geometry;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.Ar;
import net.sytes.botg.array.geometry.PolynomialSurface;
import net.sytes.botg.array.math.Mat;
import net.sytes.botg.array.math.Vec;

public class UnitTest_PolynomialSurface {

	@Test
	public void test000() {
		
		PolynomialSurface poly = new PolynomialSurface.Builder()
				.weights(new double[] { 1.0, 1.0, 1.0})
				.build();
		
		double[] x = Vec.linspace(-1.0, 1.0, 10);
		double[] y = Vec.linspace(-1.0, 1.0, 10);
		
		double[][][] XY = Mat.meshgrid(x, y);
		double[][] X = XY[0];
		double[][] Y = XY[1];
		
		poly.create(x, y);
		
		//Mat.print(poly.z());
		
	}
	
	@Test
	public void test010() {
		
		double[][] a_ij = {{1, 3}, {2, 4}};
		
		PolynomialSurface poly = new PolynomialSurface.Builder()
				.a_ij(a_ij)
				.build();
		
		System.out.println(Arrays.toString(poly.weights()));
		
	}
	
	@Test
	public void test020() {
		
		double[] weights = {1, 2, 3, 4};
		
		PolynomialSurface poly = new PolynomialSurface.Builder()
				.weights(weights)
				.build();
		
		Ar.print(poly.a_ij());
		
	}
	
	@Test
	public void test030() {
		
		double[][] a_ij = {{0, 2, 5}, {1, 3, 7}, {4, 6, 8}};
		
		PolynomialSurface poly = new PolynomialSurface.Builder()
				.a_ij(a_ij)
				.build();
		
		System.out.println(Arrays.toString(poly.weights()));
		
	}
	
}
