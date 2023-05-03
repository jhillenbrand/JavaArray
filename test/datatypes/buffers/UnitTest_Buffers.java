package datatypes.buffers;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.Ar;
import net.sytes.botg.array.math.Vec;
import net.sytes.botg.datatypes.buffers.TimedBuffer;

public class UnitTest_Buffers {
	
	@Test	
	public void test040() throws InterruptedException {
		
		TimedBuffer buf = new TimedBuffer.Builder()
				.capacity(50)
				.build();
		
		int i = 0;
		
		int s1 = 5;
		int s2 = 8;
		
		while (true) {
			
			double[] d = Vec.linspace((double) i, (double) i+1, s1);
			
			buf.push(Ar.wrap(d));
			
			System.out.println("BUF: " + Arrays.toString(buf.toArray(true)));
			
			if (buf.size() >= s2) {
				System.out.println("OUT: " + Arrays.toString(buf.toArray(false, s2)));
			}
			Thread.sleep(500);
			++i;
		}
		
	}
	
}
