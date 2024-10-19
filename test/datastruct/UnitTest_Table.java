package datastruct;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;

import net.sytes.botg.datastruct.Table;
import net.sytes.botg.datastruct.Table.CloningBehavior;
import net.sytes.botg.datastruct.Table.Comparator;

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
		
		Table t = new Table(CloningBehavior.CLONE_ON_ENTRY);
		
		Object[] a = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] b = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", a);
		t.addColumn("C2", b);
		
		a[0] = 10.0;
		
		System.out.println(t.toJson());
		
	}
	
	@Test
	public void test020() {
		
		Table t = new Table(CloningBehavior.CLONE_ON_ENTRY);
		
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
		
		t = new Table(CloningBehavior.CLONE_ON_ENTRY);
		System.out.println(t.toJson());	
				
		t = new Table("C1", Arrays.asList(c1), CloningBehavior.CLONE_ON_ENTRY);
		System.out.println(t.toJson());	
		
		t = new Table(c1);
		System.out.println(t.toJson());	
		
		t = new Table("C2", c1);
		System.out.println(t.toJson());	
		
		t = new Table(c1);
		System.out.println(t.toJson());	
		
		t = new Table("C3", Arrays.asList(c1));
		System.out.println(t.toJson());	
		
		Map<String, Object> map = new HashMap<String, Object>();
		
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
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		
		Object[] r1 = new Object[] {6.0, "F"};
		
		t.addRow(r1);
		
		System.out.println(t);
		
	}
	
	@Test
	public void test051() {
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		
		List<Object> r1 = new ArrayList<Object>();
		r1.add(6.0);
		r1.add("F");
		
		t.addRow(r1);
		
		System.out.println(t);
		
	}
	
	@Test
	public void test052() {
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		Object[] c3 = new Object[] {true, false, true, true, false};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		t.addColumn("C3", c3);
		
		Map<String, Object> map = new LinkedHashMap<String, Object>();
		map.put("C1", 6.0);
		map.put("C3", false);
		
		t.addRow(map);
		
		System.out.println(t);
		
	}
	
	@Test
	public void test060() {
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		Object[] c3 = new Object[] {true, false, true, true, false};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		t.addColumn("C3", c3);
		
		t.add("C2", "F");
		
		System.out.println(t);
		
	}
	
	@Test
	public void test061() {
		
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "C", "D", "E"};
		Object[] c3 = new Object[] {true, false, true, true, false};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		t.addColumn("C3", c3);
		
		t.add(1, "F");
		
		System.out.println(t);
		
	}
	
	@Test
	public void test070() {
		
		Table t1 = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {true, false, true, true, false};
		
		t1.addColumn("C1", c1);
		t1.addColumn("C2", c2);
				
		System.out.println(t1);
		
		Table t2 = new Table();
		Object[] c3 = new Object[] {1.0, 2.0, 3.0};
		t2.addColumn("C3", c3);
		
		System.out.println(t2);
		
		t1.addColumns(t2, true);
		
		System.out.println(t1);		
	}
	
	@Test
	public void test072() {
		
		Table t1 = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {true, false, true, true, false};
		
		t1.addColumn("C1", c1);
		t1.addColumn("C2", c2);
				
		System.out.println(t1);
		
		Table t2 = new Table();
		Object[] c3 = new Object[] {1.0, 2.0, 3.0};
		t2.addColumn("C1", c3);
		
		System.out.println(t2);
		
		t1.addRows(t2, true);
		
		System.out.println(t1);		
	}
	
	@Test
	public void test080() {
		
		Table t1 = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {true, false, true, true, false};
		
		t1.addColumn("C1", c1);
		t1.addColumn("C2", c2);
		
		System.out.println(t1);
		
		t1.set("C1", 1, 2.5);
		
		System.out.println(t1);
		
	}
	
	@Test
	public void test081() {
		
		Table t1 = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {true, false, true, true, false};
		
		t1.addColumn("C1", c1);
		t1.addColumn("C2", c2);
		
		System.out.println(t1);
		
		t1.set(1, 1, "Z");
		
		System.out.println(t1);
		
	}
	
	@Test
	public void test082() {
		
		Table t1 = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {true, false, true, true, false};
		
		t1.addColumn("C1", c1);
		t1.addColumn("C2", c2);
		
		System.out.println(t1);
		
		t1.set(1, new Object[]{6.0, "F"});
		
		System.out.println(t1);
		
	}
	
	@Test
	public void test090() {
		
		Table t1 = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {true, false, true, true, false};
		
		t1.addColumn("C1", c1);
		t1.addColumn("C2", c2);
		
		System.out.println(t1);
		
		t1.remove(1, 1);
		
		t1.remove("C1", 3);
		
		System.out.println(t1);
		
	}
	
	@Test
	public void test091() {
		
		Table t1 = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {true, false, true, true, false};
		
		t1.addColumn("C1", c1);
		t1.addColumn("C2", c2);
		
		System.out.println(t1);
		
		t1.remove("C1");
				
		System.out.println(t1);
		
	}
	
	@Test
	public void test092() {
		
		Table t1 = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 3.0, 4.0, 5.0};
		Object[] c2 = new Object[] {true, false, true, true, false};
		
		t1.addColumn("C1", c1);
		t1.addColumn("C2", c2);
		
		System.out.println(t1);
		
		t1.removeColumn(1);
				
		System.out.println(t1);
		
	}
	
	@Test
	public void test093() {
		
		Table t1 = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 1.0, 4.0, 5.0};
		Object[] c2 = new Object[] {true, false, true, true, false};
		
		t1.addColumn("C1", c1);
		t1.addColumn("C2", c2);
		
		System.out.println(t1);
		
		t1.removeDuplicates();
				
		System.out.println(t1);
		
	}
	
	@Test
	public void test094() {
		
		Table t1 = new Table();
		
		Object[] c1 = new Object[] {1.0, 2.0, 1.0, 4.0, 5.0};
		Object[] c2 = new Object[] {true, false, false, true, false};
		
		t1.addColumn("C1", c1);
		t1.addColumn("C2", c2);
		
		System.out.println(t1);
		
		t1.removeDuplicates("C1");
				
		System.out.println(t1);
		
	}
	
	@Test
	public void test095() {
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 3.0, 3.0, 3.0, 5.0};
		Object[] c2 = new Object[] {"A", "B", "B", "B", "B"};
		Object[] c3 = new Object[] {true, false, true, true, false};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		t.addColumn("C3", c3);
				
		System.out.println(t);
		
		t.removeDuplicates(new String[] {"C1", "C2"});
		
		System.out.println(t);
		
	}
	
	@Test
	public void test100() {
		
		Object obj = 1;
		
		double d = (double) ((int) obj);
		
		System.out.println(d);
	}
	
	@Test
	public void test110() {
	
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 3.3, 3.0, 3.1, 5.0};
		Object[] c2 = new Object[] {"A", "B", "B", "B", "B"};
		Object[] c3 = new Object[] {true, false, true, true, false};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		t.addColumn("C3", c3);
				
		System.out.println(t);
		
		Table t2 = t.filter(3.0, "C1", Comparator.GREATER_THAN);
		
		System.out.println(t2);
		
		Table t3 = t.filter(3.0, "C1", Comparator.EQUAL_OR_GREATER);
		
		System.out.println(t3);
		
	}	
	
	@Test
	public void test111() {
	
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 3.3, 3.0, 3.1, 5.0};
		Object[] c2 = new Object[] {"A", "B", "B", "B", "B"};
		Object[] c3 = new Object[] {true, false, true, true, false};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		t.addColumn("C3", c3);
				
		System.out.println(t);
		
		Table t2 = t.filter(3.0, "C1", Comparator.SMALLER_THAN);
		
		System.out.println(t2);
		
		Table t3 = t.filter(3.0, "C1", Comparator.EQUAL_OR_SMALLER);
		
		System.out.println(t3);
		
	}	
	
	@Test
	public void test112() {
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 3.3, 3.0, 3.1, 5.0};
		Object[] c2 = new Object[] {"A", "B", "B", "B", "C"};
		Object[] c3 = new Object[] {true, false, true, true, false};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		t.addColumn("C3", c3);
				
		System.out.println(t);
		
		Table t2 = t.filter("B", "C2", Comparator.EQUAL);
		
		System.out.println(t2);
	}
	
	@Test
	public void test113() {
		Table t = new Table();
		
		Object[] c1 = new Object[] {1.0, 3.3, 3.0, 3.1, 5.0};
		Object[] c2 = new Object[] {"A", "B", "B", "B", "C"};
		Object[] c3 = new Object[] {true, false, true, true, false};
		
		t.addColumn("C1", c1);
		t.addColumn("C2", c2);
		t.addColumn("C3", c3);
				
		System.out.println(t);
		
		Table t2 = t.filter("B", "C2", Comparator.EQUAL);
		
		System.out.println(t2);
		
		Table t3 = t2.filter(false, "C3", Comparator.EQUAL);
		
		System.out.println(t3);
	}
	
	
}
