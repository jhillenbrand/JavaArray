package geometry;

import org.junit.jupiter.api.Test;

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
	
}
