package geometry;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.Ar;
import net.sytes.botg.array.geometry.CoordSys;

public class UnitTest_CoordSys {

	@Test
	public void test000() {
		
		CoordSys cos = new CoordSys();
		
		Ar.print(cos.base());
		
		cos.rot3(45);
		
		Ar.print(cos.base());
		
		
	}
	
	@Test
	public void test001() {
		
		CoordSys cos = new CoordSys();
		
		cos.rot3(45);
		
		double[][] xyz = {{1, 0, 1},
						  {0, 1, 1},
						  {0, 0, 0}};
		
		double[][] XYZ = cos.transform(xyz);
		
		Ar.print(XYZ);
		
	}
	
	@Test
	public void test002() {
		
		CoordSys cos = new CoordSys();
		
		Ar.print(cos.base());
		
		cos.rot3(45);
		
		Ar.print(cos.base());
		
		cos.rot3(45);
		
		Ar.print(cos.base());
		
		
	}
	
}
