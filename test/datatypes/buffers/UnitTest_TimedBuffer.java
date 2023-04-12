package datatypes.buffers;

import org.junit.jupiter.api.Test;

import net.sytes.botg.datatypes.DataType;
import net.sytes.botg.datatypes.buffers.TimedBuffer;

public class UnitTest_TimedBuffer {

	@Test
	public void test000() {
		
		TimedBuffer tb = new TimedBuffer.Builder()
				.capacity(10)
				.build();
		
		for (int i = 1; i < 15; i++) {
			tb.push(System.currentTimeMillis(), i);
			System.out.println(tb);
			System.out.println("newest:" + tb.newest());
			System.out.println("oldest:" + tb.oldest());
			System.out.println(tb.toTimeSeries(true));
		}
		
	}
	
	@Test
	public void test010() {
		
		TimedBuffer tb = new TimedBuffer.Builder()
				.capacity(10)
				.build();
		
		long[] times = new long[] {1,2,3,4,5};
		Object[] elements = new Object[] {1,2,3,4,5};
		
		tb.push(times, elements);
		
		System.out.println(tb.toString());
		
	}
	
	@Test
	public void test020() {
		
		int c = 100;
		
		TimedBuffer buf = new TimedBuffer.Builder()
				.capacity(c)
				.dataType(DataType.INT)
				.build();
		
		for (int i = 0; i < c; i++) {
			
			buf.push(System.nanoTime(), i);
			
		}
		
		System.out.println(buf.toString());
		
		buf.clear(10);
		
		System.out.println(buf.toString());
		
	}
	
	@Test
	public void test030() {
		
		int c = 100;
		int j = 0;
		TimedBuffer buf = new TimedBuffer.Builder()
				.capacity(c)
				.dataType(DataType.INT)
				.build();
		
		for (int i = 0; i < c; i++) {
			
			buf.push(System.nanoTime(), j);
			++j;
		}
		
		int c1 = 60;
		
		for (int i = 0; i < c1; i++) {
			
			buf.push(System.nanoTime(), j);
			
			++j;
		}
		
		System.out.println(buf.toString());

		buf.clear(50);
		
		System.out.println(buf.toString());
	}
		
}
