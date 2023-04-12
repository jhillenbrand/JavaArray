package datatypes.graphs;

import org.junit.jupiter.api.Test;

import net.sytes.botg.datatypes.graphs.DirectedGraphNode;

public class UnitTest_DirectedGraphNode {

	@Test
	public void test000() {
		DirectedGraphNode n1 = new DirectedGraphNode("p1");
		
		DirectedGraphNode n2 = new DirectedGraphNode("p2");
		
		n1.add(n2);
		
		DirectedGraphNode n3 = new DirectedGraphNode("p3");
		
		n2.add(n3);
		
		DirectedGraphNode n4 = new DirectedGraphNode("p4");
		
		n2.add(n4);
		
		n4.add(n1);
		
		DirectedGraphNode root = n4.findStartNode(null);
				
		System.out.println(root.toString());
	}
	
	@Test
	public void test010() {
		
		DirectedGraphNode n1 = new DirectedGraphNode("N1");
		DirectedGraphNode n2 = new DirectedGraphNode("N2");
		DirectedGraphNode n3 = new DirectedGraphNode("N3");
		DirectedGraphNode n4 = new DirectedGraphNode("N4");
		DirectedGraphNode n5 = new DirectedGraphNode("N5");
		DirectedGraphNode n6 = new DirectedGraphNode("N6");
		
		n1.add(n2);
		n2.add(n3);
		n3.add(n4);
		n4.add(n5);
		n5.add(n6);
		n6.add(n2); 
		
		DirectedGraphNode sn = n6.findStartNode(null);
		System.out.println(sn.toString());
		
	}
	
}
