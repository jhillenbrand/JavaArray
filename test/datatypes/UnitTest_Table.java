package datatypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
	
	@Test
	public void test010() {
		
		Table t = new Table(true);
		
		Object[] a = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] b = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", a);
		t.addColumn("C2", b);
		
		a[0] = 10.0;
		
		System.out.println(t.toJson());
		
	}
	
	@Test
	public void test020() {
		
		Table t = new Table(false);
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		
		Object[] r1 = new Object[] {6.0, "F"};
		
		t.addRow(r1);
		
		r1[0] = 10.0;
		
		System.out.println(t.toJson());		
		
	}

	@Test
	public void test030() {
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		Table t = new Table();
		System.out.println(t.toJson());	
		
		t = new Table(true);
		System.out.println(t.toJson());	
				
		t = new Table("C1", Arrays.asList(c1), true);
		System.out.println(t.toJson());	
		
		t = new Table(c1);
		System.out.println(t.toJson());	
		
		t = new Table("C2", c1);
		System.out.println(t.toJson());	
		
		t = new Table(Arrays.asList(c1));
		System.out.println(t.toJson());	
		
		t = new Table("C3", Arrays.asList(c1));
		System.out.println(t.toJson());	
		
		Map<String, List<Object>> map = new HashMap<String, List<Object>>();
		
		map.put("C4", Arrays.asList(c1));
		map.put("C5", Arrays.asList(c2));
		
		t = new Table(map);
		System.out.println(t.toJson());	
		
		Object[][] data = new Object[2][];
		data[0] = c1;
		data[1] = c2;
		
		t = new Table(data);
		System.out.println(t.toJson());	
		
		String[] headers = new String[] {"C6", "C7"};
		t = new Table(headers, data);
		System.out.println(t.toJson());	
		
		List<String> headers2 = new ArrayList<String>();
		headers2.add("C8");
		headers2.add("C9");
		t = new Table(headers2, data);
		System.out.println(t.toJson());
		
		t = new Table(headers);
		System.out.println(t.toJson());
	}
		
	@Test
	public void test040() {
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		
		System.out.println(t.get(0, 0));
		
		System.out.println(t.get(1, 1));
		
		//System.out.println(t.get(1, 2));
		
	}
	
	@Test
	public void test041() {
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		
		System.out.println(t.getRow(1));
		
	}
	
	@Test
	public void test042() {
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);

		System.out.println(t.getColumn("C1"));
		System.out.println(t.getColumn(1));
		
		List<Object> list = t.getColumn("C1");
		
		list.set(0, 10.0);
		
		System.out.println(list);
		
		System.out.println(t);		
		
	}
	
	@Test
	public void test043() {
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);

		System.out.println(Arrays.toString(t.getColumnNames()));
		
	}
	
	@Test
	public void test044() {
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);

		System.out.println(t.getColumnName(2));
	}
	
	@Test
	public void test045() {
		
		Table t = new Table();
		
		//Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		//Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1");

		System.out.println(t);
		
	}
	
	@Test
	public void test046() {
		
		Table t = new Table();
		
		//Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		//Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1");
		t.addColumn("C2");

		System.out.println(t);
		
	}
	
	@Test
	public void test047() {
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		//Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", c1);
		
		t.addColumn("C2");

		System.out.println(t);
		
	}
	
	@Test
	public void test050() {
		
	}
	
}
