package math;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.math.Mat;
import net.sytes.botg.array.math.Vec;

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
	
		Mat.print(points);
		
		double[][] D_ij = Mat.distanceMatrix(points);
		
		Mat.print(D_ij);
		
	}
	
	@Test
	public void test010() {
		
		double[][] m1 = new double[2][3];
		double[][] m2 = new double[1][3];
		
		m1[0] = new double[] {0.0, 0.1, 0.2};
		m1[1] = new double[] {1.0, 1.1, 1.2};
		
		m2[0] = new double[] {2.0, 2.1, 2.2};
		
		double[][] m3 = Mat.appendRows(m1, m2);
		
		Mat.print(m3);
	}
	
	@Test
	public void test020() {
		
		double[][] m1 = new double[2][3];
		double[][] m2 = new double[2][4];
		
		m1[0] = new double[] {0.0, 0.1, 0.2};
		m1[1] = new double[] {1.0, 1.1, 1.2};
		
		m2[0] = new double[] {2.0, 2.1, 2.2, 2.3};
		m2[1] = new double[] {3.0, 3.1, 3.2, 3.3};
		
		double[][] m3 = Mat.appendRows(m1, m2);
		
		Mat.print(m3);
	}
	
	@Test
	public void test030() {
		
		double[][] m1 = new double[2][3];
		
		m1[0] = new double[] {0.0, 0.1, 0.2};
		m1[1] = new double[] {1.0, 1.1, 1.2};
		
		double[][] m2 = Mat.transpose(m1);
		
		Mat.print(m2);
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
		
		double[][] X = Mat.zeros(2, 2);
		
		X[0][0] = 2;
		X[1][0] = 6;
		X[0][1] = 1;
		X[1][1] = 4;
		
		double[][] X_inv = Mat.inverse(X);
		
		Mat.print(X);
		Mat.print(X_inv);
		
	}
	
	@Test
	public void test060() {
		
		double[][] X = Mat.incrementMat(3, 3);
				
		double[][] Y = Mat.subByElimination(X, 1, 1);
		
		Mat.print(X);
		Mat.print(Y);
		
	}
	
	@Test
	public void test070() {
	
		double[][] X  = {{0, 1, 2},{3, 2, 1},{1, 1, 0}};
		
		double d = Mat.det(X);
		
		assertEquals(d, 3);
		
	}
	
	@Test
	public void test080() {
		
		double[][] I = Mat.unitMatrix(3);
		
		double[][] I_inv = Mat.inverse(I);
		
		Mat.print(I);
		Mat.print(I_inv);
		
	}
	
	@Test
	public void test081() {
		
		double[][] X = {{1, 2}, {2, 3}};
		
		double[][] X_inv = Mat.inverse(X);
		
		Mat.print(X);
		Mat.print(X_inv);
		
	}
	
	@Test
	public void test082() {
		
		double[][] X = {{1, 2, 0}, {2, 4, 1}, {2, 1, 0}};
		
		double[][] X_inv = Mat.inverse(X);
		
		Mat.print(X);
		Mat.print(X_inv);
		
 	}
	
	@Test
	public void test083() {
		double[][] X = {{1, 2, 0}, {2, 4, 1}, {2, 1, 0}};
		
		double[][] X_inv = Mat.inverse(X);
		
		Mat.print(X);
		Mat.print(X_inv);
		
		double[][] X_inv_appr = Mat.inverseNumerical(X, 100);
		Mat.print(X_inv_appr);
	}
	
	@Test
	public void test084() {
		double[][] X = {{1, 2}, {2, 1}};
		
		double[][] X_inv = Mat.inverse(X);
		
		Mat.print(X);
		Mat.print(X_inv);
		
		double[][] X_inv_appr = Mat.inverseNumerical(X, 100);
		Mat.print(X_inv_appr);
	}
	
	@Test
	public void test090() {
		double[][] X = {{0, 0},{1, 0}};
		double[][] Y = {{0, 1}, {0, 0}};
		
		double[][] Z = Mat.product(X, Y);
		double[][] Z2 = Mat.product(X, Y);
		
		Mat.print(Z);
		Mat.print(Z2);
	}
	
	@Test
	public void test091() {
		double[][] X = {{1, 0, 1}, {2, 1, 1}, {0, 1, 1}, {1, 1, 2}};
		double[][] Y = {{1, 2, 1}, {2, 3, 1}, {4, 2, 2}};
		
		double[][] Z = Mat.product(X, Y);
		double[][] Z2 = Mat.product2(X, Y);
		double[][] Z3 = Mat.product3(X, Y);
		
		Mat.print(Z);
		Mat.print(Z2);
		Mat.print(Z3);
		
	}	
	
	@Test
	public void test092() {
		double[][] X = {{1, 1}, {1, 1}};
		
		double[][] X_2 = Mat.square(X);
		double[][] X_3 = Mat.power(X, 3);
		
		Mat.print(X);
		Mat.print(X_2);
		Mat.print(X_3);
		
	}
	
	@Test
	public void test093() {
		double[][] X = Mat.ones(3, 3);
		
		double[][] X_2 = Mat.square(X);
		double[][] X_3 = Mat.power(X, 3);
		
		Mat.print(X);
		Mat.print(X_2);
		Mat.print(X_3);
		
	}
	
	@Test
	public void test094() {
		
		long st = 0;
		long et = 0;
		int n = 500;
		int r = 10;

		double[][] A = Mat.ones(n, n);
		double[][] B = Mat.ones(n, n);
		
		st = System.nanoTime();
		
		for (int i = 0; i < r; i++) {
		
			
			double[][] C = Mat.product2(A, B);
			
			//Mat.print(C);
			
		}
		
		et = System.nanoTime();
		
		System.out.println("NAIVE n=" + n);
		System.out.println("elapsed time [ns]: " + (et - st) + ", per iteration [ns]: " + (double) (et - st) /  (double) r);
		System.out.println("");
		
		st = System.nanoTime();
		
		for (int i = 0; i < r; i++) {
		
			
			double[][] C = Mat.product3(A, B);
			
			//Mat.print(C);
			
		}
		
		et = System.nanoTime();
		
		System.out.println("ROW-ORIENTED n=" + n);
		System.out.println("elapsed time [ns]: " + (et - st) + ", per iteration [ns]: " + (double) (et - st) /  (double) r);
		System.out.println("");
		
		st = System.nanoTime();
		
		for (int i = 0; i < r; i++) {
		
			
			double[][] C = Mat.product(A, B);
			
			//Mat.print(C);
			
		}
		
		et = System.nanoTime();
		
		System.out.println("JAMA n=" + n);
		System.out.println("elapsed time [ns]: " + (et - st) + ", per iteration [ns]: " + (double) (et - st) /  (double) r);
		System.out.println("");
	}
	
	@Test
	public void test100() {
		double[][] A = {{2, 1}, {1, 2}};
		double[] y = {1, 1};
		
		double[] x = Mat.linSolveInverse(A, y);
		double[] x2 = Mat.linSolveLU(A, y);
		
		Vec.print(x);
		Vec.print(x2);
	}
	
	@Test
	public void test101() {
		double[][] A = {{0, 1, 3, 0, 1}, {0, 2, 0, 0, 4}, {3, 0, 2, 0, 0}, {3, 1, 0, 1, 0}, {0, 1, 1, 0, 1}};
		double[] y = {1, 1, 0, 4, 2};
		
		double[] x = Mat.linSolveInverse(A, y);
		double[] x2 = Mat.linSolveLU(A, y);
		double[] x3 = Mat.linSolveGaussian(A, y);
		Vec.print(x);
		Vec.print(x2);
		Vec.print(x3);
	}
	
	@Test
	public void test102() {
		double[][] A = {{0, 1, 3, 0}, {2, 0, 0, 4}, {3, 0, 2, 0}, {3, 1, 0, 1}};
		double[] y = {1, 1, 0, 4};
		
		double[] x = Mat.linSolveInverse(A, y);		
		double[] x2 = Mat.linSolveLU(A, y);
		double[] x3 = Mat.linSolveGaussian(A, y);
		Vec.print(x);
		Vec.print(x2);
		Vec.print(x3);
	}
		
	@Test
	public void test103() {
		
		long st = 0;
		long et = 0;
		int n = 20;
		
		for (int i = 2; i < n; i++) {
		
			double[][] A = Mat.ones(i, i);
			double[] y = Vec.ones(i);
			
			st = System.currentTimeMillis();
			double[] x = Mat.linSolveInverse(A, y);
			et = System.currentTimeMillis();
			
			//Mat.print(x);
			System.out.println("n=" + i);
			System.out.println("elapsed time [ms]: " + (et - st));
			System.out.println("");
		}
	}
	
	@Test
	public void test104() {
		
		long st = 0;
		long et = 0;
		int n = 20;
		
		for (int i = 2; i < n; i++) {
		
			double[][] A = Mat.unitMatrix(i);
			double[] y = Vec.ones(i);
			
			st = System.currentTimeMillis();
			double[] x = Mat.linSolveInverse(A, y);
			et = System.currentTimeMillis();
			
			Vec.print(x);
			System.out.println("n=" + i);
			System.out.println("elapsed time [ms]: " + (et - st));
			System.out.println("");
		}
	}
	
	@Test
	public void test105() {
		
		long st = 0;
		long et = 0;
		int n = 20;
		
		for (int i = 2; i < n; i++) {
		
			double[][] A = Mat.unitMatrix(i);
			double[] y = Vec.ones(i);
			
			st = System.currentTimeMillis();
			double[] x = Mat.linSolveGaussian(A, y);
			et = System.currentTimeMillis();
			
			Vec.print(x);
			System.out.println("n=" + i);
			System.out.println("elapsed time [ms]: " + (et - st));
			System.out.println("");
		}
	}
	
	@Test
	public void test106() {
		double[][] A = {{3.0, 2.0, -4.0},{2.0, 3.0, 3.0},{5.0, -3.0, 1.0}};
		double[] y = {3.0, 15.0, 14.0};
		
		double[] x = Mat.linSolveLU(A, y);
		
		Vec.print(x);
		
	}
	
	@Test
	public void test150() {
		double[][] X = Mat.incrementMat(3,  3);
		Mat.print(X);
	}
	
	@Test
	public void test151() {
		double[][] X = Mat.incrementMat(2,  5);
		Mat.print(X);
	}
	
	@Test
	public void test152() {
		double[][] X = Mat.incrementMat(5,  3);
		Mat.print(X);
	}

	@Test
	public void test134() {
		double[][] nans = Mat.nan(10, 2);
		
		System.out.println(Arrays.deepToString(nans));
		
	}
	
	
	@Test
	public void test130() {
		
		int n = 100;
		int r = 20_000;
		
		double[][] X = Mat.unitMatrix(n);
		double[][] Y = null;
		
		long st = System.nanoTime();
		for (int i = 0; i < r; i++) {
			
			Y = Mat.copy2(X);
			
		}
		
		long et = System.nanoTime();
		double el = et - st;
		double sp = el / n; 
		
		System.out.println("copy2 -> Sampling Period per Element [ns]: " + sp);
		
		st = System.nanoTime();
		for (int i = 0; i < r; i++) {
			
			Y = Mat.copy(X);
			
		}
		
		et = System.nanoTime();
		el = et - st;
		sp = el / n; 
		
		System.out.println("copy -> Sampling Period per Element [ns]: " + sp);
		
		System.out.println("X = [" + X[0].length + ", " + X.length + "], Y = [" + Y[0].length + ", " + Y.length + "]");
		
	}

	@Test
	public void test131() {
		
		double[][] X = Mat.unitMatrix(3);
		
		double[][] Y = Mat.copy(X);
		
		Y[0][0] = 2.0;
		
		Mat.print(X);
		Mat.print(Y);
		
	}
	
	
	@Test
	public void test074() {
		
		double[] x = Vec.linspace(-1.0, 1.0, 10);
		double[] y = Vec.linspace(-1.0, 1.0, 5);
		
		double[][][] XY = Mat.meshgrid(x, y);
				
		Mat.print(XY[0]);
		Mat.print(XY[1]);
		
	}
}
