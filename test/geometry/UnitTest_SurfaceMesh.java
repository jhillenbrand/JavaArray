package geometry;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.Ar;
import net.sytes.botg.array.geometry.Paraboloid;
import net.sytes.botg.array.math.Mat;
import net.sytes.botg.array.math.Vec;

public class UnitTest_SurfaceMesh {

	@Test
	public void test000() {
		
		Paraboloid p = new Paraboloid.Builder()
				.a(1)
				.b(1)
				.elliptic(true)
				.build();
		
		double[] x = Vec.linspace(-1.0, 1.0, 100);
		double[] y = Vec.linspace(-1.0, 1.0, 100);
		
		p.create(x, y);
		
		Ar.print(p.x());
		Ar.print(p.y());
		
		Ar.print(p.z());
		
		
	}
	
}
