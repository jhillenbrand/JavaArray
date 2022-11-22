package math;

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
		
		m2[0] = new double[] {2.0, 2.1, 2.2, 2.3};
		
		double[][] m3 = Mat.appendByColumn(m1, m2);
		
		ArUtils.print(m3);
	}
	
}
