package datatypes.buffers;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import net.sytes.botg.datatypes.DataType;
import net.sytes.botg.datatypes.buffers.DuplicateTimedBuffer;

public class UnitTest_DuplicateTimedBuffer {

	@Test
	public void test000() {
		 
		DuplicateTimedBuffer dtb = new DuplicateTimedBuffer.Builder()
				.capacity(3)
				.dataType(DataType.DOUBLE)
				.id("B1")
				.duplicateIds(Arrays.asList(new String[]{"B2", "B3"}))
				.build();
				
		for (int i = 0; i < 100; i++) {
			
			dtb.push((double) i);
						
			System.out.println(dtb.toString());			
			System.out.println(dtb.getDuplicates().toString());
			
  			System.out.println(i);
			
		}
		
	}
	
}
