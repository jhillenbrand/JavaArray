package math;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.jupiter.api.Test;

import net.sytes.botg.array.ArUtils;
import net.sytes.botg.array.math.Mat;

public class UnitTest_Mat {

	@Test
	public void test000() {
		
		double[][] points = new double[5][2];
		
		points[0][0] = 1.0;
		points[0][1] = 1.0;
		points[1][0] = 2.0;
		points[1][1] = 1.0;
		points[2][0] = 4.0;
		points[2][1] = 3.0;
		points[3][0] = 3.0;
		points[3][1] = 2.0;
		points[4][0] = 3.0;
		points[4][1] = 3.0;
	
		ArUtils.print(points);
		
		double[][] D_ij = Mat.distanceMatrix(points);
		
		ArUtils.print(D_ij);
		
	}
	
	@Test
	public void test010() {
		
		double[][] m1 = new double[2][3];
		double[][] m2 = new double[1][3];
		
		m1[0] = new double[] {0.0, 0.1, 0.2};
		m1[1] = new double[] {1.0, 1.1, 1.2};
		
		m2[0] = new double[] {2.0, 2.1, 2.2};
		
		double[][] m3 = Mat.appendByColumns(m1, m2);
		
		ArUtils.print(m3);
	}
	
	@Test
	public void test020() {
		
		double[][] m1 = new double[2][3];
		double[][] m2 = new double[2][4];
		
		m1[0] = new double[] {0.0, 0.1, 0.2};
		m1[1] = new double[] {1.0, 1.1, 1.2};
		
		m2[0] = new double[] {2.0, 2.1, 2.2, 2.3};
		m2[1] = new double[] {3.0, 3.1, 3.2, 3.3};
		
		double[][] m3 = Mat.appendByRows(m1, m2);
		
		ArUtils.print(m3);
	}
	
	@Test
	public void test030() {
		
		double[][] m1 = new double[2][3];
		
		m1[0] = new double[] {0.0, 0.1, 0.2};
		m1[1] = new double[] {1.0, 1.1, 1.2};
		
		double[][] m2 = Mat.transpose(m1);
		
		ArUtils.print(m2);
	}
	
	@Test
	public void test040() {
		
		double[][] M = new double[2][2];
		
		M[0][0] = 1.0;
		M[0][1] = 1.0;
		M[1][0] = 1.0;
		M[1][1] = 1.0;
		
		double[][] N = M.clone();
		
		N[0][1] = 2.0;
		
		System.out.println(Arrays.deepToString(M));
		
	}
	
	@Test
	public void test050() {
		
		double[][] X = ArUtils.zeros(2, 2);
		
		X[0][0] = 2;
		X[1][0] = 6;
		X[0][1] = 1;
		X[1][1] = 4;
		
		double[][] X_inv = Mat.inverse(X);
		
		ArUtils.print(X);
		ArUtils.print(X_inv);
		
	}
	
	@Test
	public void test060() {
		
		double[][] X = ArUtils.incrementMat(3, 3);
				
		double[][] Y = Mat.sub(X, 1, 1);
		
		ArUtils.print(X);
		ArUtils.print(Y);
		
	}
	
	@Test
	public void test070() {
	
		double[][] X  = {{0, 1, 2},{3, 2, 1},{1, 1, 0}};
		
		double d = Mat.det(X);
		
		assertEquals(d, 3);
		
	}
	
	@Test
	public void test080() {
		
		double[][] I = ArUtils.unitMatrix(3);
		
		double[][] I_inv = Mat.inverse(I);
		
		ArUtils.print(I);
		ArUtils.print(I_inv);
		
	}
	
	@Test
	public void test081() {
		
		double[][] X = {{1, 2}, {2, 3}};
		
		double[][] X_inv = Mat.inverse(X);
		
		ArUtils.print(X);
		ArUtils.print(X_inv);
		
	}
	
	@Test
	public void test082() {
		
		double[][] X = {{1, 2, 0}, {2, 4, 1}, {2, 1, 0}};
		
		double[][] X_inv = Mat.inverse(X);
		
		ArUtils.print(X);
		ArUtils.print(X_inv);
		
 	}
	
	@Test
	public void test090() {
		double[][] X = {{0, 0},{1, 0}};
		double[][] Y = {{0, 1}, {0, 0}};
		
		double[][] Z = Mat.prod(X, Y);
		
		ArUtils.print(Z);
		
	}
	
	@Test
	public void test091() {
		double[][] X = {{1, 0, 1}, {2, 1, 1}, {0, 1, 1}, {1, 1, 2}};
		double[][] Y = {{1, 2, 1}, {2, 3, 1}, {4, 2, 2}};
		
		double[][] Z = Mat.prod(X, Y);
		
		ArUtils.print(Z);
		
	}
	
	@Test
	public void test100() {
		double[][] A = {{2, 1}, {0, 2}};
		double[] y = {1, 1};
		
		double[] x = Mat.linSolve(A, y);
		
		ArUtils.print(x);
	}
	
}
