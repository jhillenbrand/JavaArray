package testing.MathArray;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import net.sytes.botg.array.MathArray;

public class UnitTest_MathArray {

	@Test
	public void testClosestExponentForBase2() {
		
		assertEquals(MathArray.closestExponentForBase2(1024), 10);
		
		assertEquals(MathArray.closestExponentForBase2(2047), 10);
		
	}
	
}
