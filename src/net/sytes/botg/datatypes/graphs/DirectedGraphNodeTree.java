package net.sytes.botg.datatypes.graphs;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

public class DirectedGraphNodeTree<K, E> {

	private DirectedGraphNode<E> startingNode = null;
	
	private Map<K, DirectedGraphNode<E>> treeNodes = new HashMap<K, DirectedGraphNode<E>>();
	
	public DirectedGraphNodeTree() {
		
	}
	
	public DirectedGraphNode<E> get(K key) {
		return this.treeNodes.get(key);
	}
	
	public void add(K key, DirectedGraphNode<E> parentNode, DirectedGraphNode<E> childNode) {
		if (parentNode != null) {
			if (this.treeNodes.containsValue(parentNode)) {
				parentNode.add(childNode);
			} else {
				throw new IllegalArgumentException(parentNode.toString() + " must already be part of " + this.getClass().getSimpleName());
			}
		} else {
			// store starting node
			this.startingNode = childNode;
		}
		this.treeNodes.put(key, childNode);		
	}
	
	public boolean containsKey(K key) {
		return this.treeNodes.containsKey(key);
	}
	
	public int size() {
		return this.treeNodes.size();
	}
	
	public boolean isEmpty() {
		return this.treeNodes.isEmpty();
	}
	
	public Set<K> keySet(){
		return this.treeNodes.keySet();
	}
	
	public Set<Entry<K, DirectedGraphNode<E>>> entrySet() {
		return this.treeNodes.entrySet();
	}
	
	public List<E> values(){
		return (List<E>) this.treeNodes.values();
	}
	
	public DirectedGraphNode<E> startingNode(){
		return this.startingNode;
	}
	
}
