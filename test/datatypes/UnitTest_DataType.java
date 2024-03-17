package datatypes;

import static org.junit.Assert.assertEquals;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import net.sytes.botg.array.Ar;
import net.sytes.botg.datatypes.DataType;

public class UnitTest_DataType {

	@Test
	public void test000() {		
		assertEquals(1.0, DataType.cast(true, DataType.BOOLEAN, DataType.DOUBLE));		
	}
	
	@Test
	public void test001() {		
		assertEquals(0.0, DataType.cast(false, DataType.BOOLEAN, DataType.DOUBLE));		
	}
	
	@Test
	public void test002() {		
		assertEquals(1L, DataType.cast(true, DataType.BOOLEAN, DataType.LONG));		
	}
	
	@Test
	public void test003() {		
		assertEquals(0L, DataType.cast(false, DataType.BOOLEAN, DataType.LONG));		
	}
	
	@Test
	public void test004() {		
		assertEquals((byte) 0, DataType.cast(false, DataType.BOOLEAN, DataType.BYTE));		
	}
	
	@Test
	public void test005() {		
		assertEquals("true", DataType.cast(true, DataType.BOOLEAN, DataType.STRING));		
	}
	
	@Test
	public void test006() {		
		assertEquals(100.0, DataType.cast(100L, DataType.LONG, DataType.DOUBLE));		
	}
	
	@Test
	public void test007() {		
		assertEquals(100.0f, DataType.cast(100L, DataType.LONG, DataType.FLOAT));		
	}
	
	@Test
	public void test008() {		
		assertEquals(100L, DataType.cast(100L, DataType.LONG, DataType.LONG));		
	}
	
	@Test
	public void test009() {		
		assertEquals(100, DataType.cast(100L, DataType.LONG, DataType.INT));		
	}
	
	@Test
	public void test010() {
		Double d = 100.0;
		assertEquals(100L, DataType.cast(d, DataType.DOUBLE, DataType.LONG));		
	}
	
	@Test
	public void test011() {
		short s = 11;
		assertEquals(11.0, DataType.cast(s, DataType.SHORT, DataType.DOUBLE));		
	}
	
	@Test
	public void test012() {
		Short s = 11;
		assertEquals(11L, DataType.cast(s, DataType.SHORT, DataType.LONG));		
	}
	
	@Test
	public void test013() {
		Short s = 11;
		assertEquals("11", DataType.cast(s, DataType.SHORT, DataType.STRING));		
	}
	
	@Test
	public void test014() {
		long l = 11L;
		assertEquals((short) 11, DataType.cast(l, DataType.LONG, DataType.SHORT));		
	}
	
	@Test
	public void test015() {
		assertEquals(DataType.of("true"), DataType.BOOLEAN);
	}
	
	@Test
	public void test016() {
		assertEquals(DataType.of("ja"), DataType.BOOLEAN);
	}
	
	@Test
	public void test017() {
		assertEquals(DataType.of("1.05"), DataType.DOUBLE);
	}
	
	@Test
	public void test018() {
		assertEquals(DataType.of("1.05a"), DataType.DOUBLE);
	}
	
	@Test
	public void test019() {
		assertEquals(DataType.of("1.0"), DataType.INT);
	}
	
	@Test
	public void test020() {
		assertEquals(DataType.of("124236523452342134"), DataType.LONG);
	}
	
	@Test
	public void test021() {
		assertEquals(DataType.of("12342134"), DataType.INT);
	}
	
	@Test
	public void test022() {
		assertEquals(DataType.of("123.12d"), DataType.DOUBLE);
	}
	
	@Test
	public void test023() {
		
		Map<String, Double> map = new HashMap<String, Double>();
		map.put("1", 1.0);
		map.put("2", 10.1);
		
		Map<String, Integer> map2 = DataType.castMapValueType(map, DataType.INT);
		System.out.println(map2);
	}
	
	@Test
	public void test030() {
		
		System.out.println(DataType.cast("99.9", DataType.SHORT));
		System.out.println(DataType.cast("100.0", DataType.INT));
		System.out.println(DataType.cast("100.0", DataType.LONG));
		
	}
	
	@Test
	public void test031() {
		
		String s = "";
		double d = 1.1;
		int in = 0;
		
		long t1 = System.nanoTime();
		int n = 1_000_000;
		for (int i = 0; i < n; i++) {
			s = "" + d * i;
			in = (int) DataType.cast(s, DataType.INT);
		}
		
		long t2 = System.nanoTime();
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		
	}
	
}
