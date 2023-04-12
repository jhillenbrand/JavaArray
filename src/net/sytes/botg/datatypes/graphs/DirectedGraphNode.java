package net.sytes.botg.datatypes.graphs;

import java.util.ArrayList;
import java.util.List;
import java.util.UUID;

public class DirectedGraphNode<T> {

	private List<DirectedGraphNode<T>> parents;
	private List<DirectedGraphNode<T>> children;
	private T value;
	private String id;
	
	private int level;

	public DirectedGraphNode(String id, T value) {
		this.id = id;
		this.parents = new ArrayList<DirectedGraphNode<T>>();
		this.children = new ArrayList<DirectedGraphNode<T>>();
		this.value = value;
		this.level = 0;
	}
	
	public DirectedGraphNode(T value) {
		this(UUID.randomUUID().toString(), value);
	}
	
	public DirectedGraphNode(String id) {
		this(id, null);
	}
	
	public DirectedGraphNode() {
		this(UUID.randomUUID().toString(), null);
	}
	
	public void add(DirectedGraphNode<T> node) {
		if (node == null) {
			throw new IllegalArgumentException(DirectedGraphNode.class.getSimpleName() + " must not be null");
		}
		if (node.level == 0) {
			node.level = this.level + 1; 
		}
		this.children.add(node);
		node.getParents().add(this);
	}

	private boolean hasParents() {
		if (this.parents.size() > 0) {
			return true;
		} else {
			return false;
		}
	}
	 
	private boolean hasChildren() {
		return this.children.size() > 0;
	}
	
	public DirectedGraphNode<T> findStartNode(DirectedGraphNode<T> foundNode){
		if (this.level == 0) {
			foundNode = this;
		}
		if (foundNode == null) {
			for (DirectedGraphNode<T> parentNode : this.parents) {
				foundNode = parentNode.findStartNode(foundNode);
				if (foundNode != null) {
					return foundNode;
				}
			}
			return foundNode;
		} else {
			return foundNode;
		}
	}
	
	public T getValue() {
		return this.value;
	}

	public List<DirectedGraphNode<T>> getParents() {
		return this.parents;
	}

	public List<DirectedGraphNode<T>> getChildren() {
		return this.children;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(this.getClass().getSimpleName())
			.append(" Id: ")
			.append(this.id)
			.append(", Level: ")
			.append(this.level)
			.append("\n\t");
		
		if (this.value != null) {
			sb.append(this.value.toString());
		}			
		return sb.toString();		
	}
}
