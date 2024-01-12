package math;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.math.Mat;
import net.sytes.botg.array.math.Parallel;

public class UnitTest_ParallelComputing {


	@Test
	public void test000() throws Exception {
		
		long st = 0;
		long et = 0;
		int n = 2048;
		int r = 2;

		double[][] A = Mat.rand(n, n);
		double[][] B = Mat.rand(n, n);
		
		st = System.nanoTime();
		
		for (int i = 0; i < r; i++) {
		
			
			double[][] C = Mat.product2(A, B);
			
			//Ar.print(C);
			
		}
		
		et = System.nanoTime();
		
		System.out.println("NAIVE i-> j -> k, n=" + n);
		System.out.println("elapsed time [ns]: " + (et - st) + ", per iteration [ns]: " + (double) (et - st) /  (double) r);
		System.out.println("");
		
		st = System.nanoTime();
		
		for (int i = 0; i < r; i++) {
		
			
			double[][] C = Mat.product3(A, B);
			
			//Ar.print(C);
			
		}
		
		et = System.nanoTime();
		
		System.out.println("ROW-ORIENTED n=" + n);
		System.out.println("elapsed time [ns]: " + (et - st) + ", per iteration [ns]: " + (double) (et - st) /  (double) r);
		System.out.println("");
		
		st = System.nanoTime();
		
		for (int i = 0; i < r; i++) {
		
			
			double[][] C = Mat.product(A, B);
			
			//Ar.print(C);
			
		}
		
		et = System.nanoTime();
		
		System.out.println("JAMA n=" + n);
		System.out.println("elapsed time [ns]: " + (et - st) + ", per iteration [ns]: " + (double) (et - st) /  (double) r);
		System.out.println("");
		
		st = System.nanoTime();
		
		for (int i = 0; i < r; i++) {
		
			
			double[][] C = Mat.product4(A, B);
			
			//Ar.print(C);
			
		}
		
		et = System.nanoTime();
		
		System.out.println("NAIVE i -> k -> j,  n=" + n);
		System.out.println("elapsed time [ns]: " + (et - st) + ", per iteration [ns]: " + (double) (et - st) /  (double) r);
		System.out.println("");
						
		st = System.nanoTime();
		
		for (int i = 0; i < r; i++) {
		
			
			double[][] C = Parallel.product2(A, B);
			
			//Ar.print(C);
			
		}
		
		et = System.nanoTime();
		
		System.out.println("Parallel NAIVE i -> j -> k,  n=" + n);
		System.out.println("elapsed time [ns]: " + (et - st) + ", per iteration [ns]: " + (double) (et - st) /  (double) r);
		System.out.println("");
		
		
		st = System.nanoTime();
		
		for (int i = 0; i < r; i++) {
		
			
			double[][] C = Parallel.product(A, B);
			
			//Ar.print(C);
			
		}
		
		et = System.nanoTime();
		
		System.out.println("Parallel Tiled i -> j -> k,  n=" + n);
		System.out.println("elapsed time [ns]: " + (et - st) + ", per iteration [ns]: " + (double) (et - st) /  (double) r);
		System.out.println("");
		
	}
	
}
