package datatypes;

import org.junit.jupiter.api.Test;

import net.sytes.botg.datatypes.Table;

public class UnitTest_Table {

	@Test
	public void test000() {
		
		Table t = new Table();
		
		Object[] a = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] b = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", a);
		t.addColumn("C2", b);
		
		System.out.println(t.toJson());
		
	}
	
}
