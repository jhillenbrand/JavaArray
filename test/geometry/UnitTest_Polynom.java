package geometry;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.Ar;
import net.sytes.botg.array.geometry.Polynom;
import net.sytes.botg.array.math.Vec;

public class UnitTest_Polynom {

	@Test
	public void test000() {
		Polynom poly = new Polynom.Builder()
				.weights(new double[] { 1.0, 1.0, 1.0})
				.build();
		
		double[] t = Vec.linspace(0.0, 4.0, 50);
		
		poly.create(t);
		
		Ar.print(poly.y());
		
				
	}
	
}
