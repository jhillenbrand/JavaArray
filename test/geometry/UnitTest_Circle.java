package geometry;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.ArUtils;
import net.sytes.botg.array.geometry.Circle;

public class UnitTest_Circle {

	@Test
	public void test000() {
		
		Circle c = new Circle.Builder()
				.D(20.0)
				.x_0(10.0)
				.y_0(5.0)
				.build();
		
		double[] t = ArUtils.linspace(0.0, 1.0, 100);
		
		c.create(t);
		
		double[] x = c.x();
		double[] y = c.y();
		
		System.out.println(Arrays.toString(x));		
		System.out.println(Arrays.toString(y));
	}
	
}
