package net.sytes.botg.datatypes.graphs;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * code is based on <a href="https://github.com/gt4dev/yet-another-tree-structure">https://github.com/gt4dev/yet-another-tree-structure</a>
 * @author hillenbrand
 *
 * @param <T>
 */
public class TreeNode<T> implements Iterable<TreeNode<T>> {

	private T value;
    private TreeNode<T> parent;
    private List<TreeNode<T>> children;

    public TreeNode(T data) {
        this.value = data;
        this.children = new LinkedList<TreeNode<T>>();
    }

    public TreeNode<T> addChild(T child) {
        TreeNode<T> childNode = new TreeNode<T>(child);
        childNode.parent = this;
        this.children.add(childNode);
        return childNode;
    }
    
    public boolean isRoot() {
		return this.parent == null;
	}

	public boolean isLeaf() {
		return this.children.size() == 0;
	}
	
	@Override
	public Iterator<TreeNode<T>> iterator() {
		TreeNodeIter<T> iter = new TreeNodeIter<T>(this);
		return iter;
	}

	public T getValue() {
		return this.value;
	}

	private enum ProcessStages {
		PROCESS_PARENT, PROCESS_CHILD_CURRENT_NODE, PROCESS_CHILD_SUB_NODE
	}
	
	private class TreeNodeIter<T> implements Iterator<TreeNode<T>> {
		
		private TreeNode<T> treeNode;
		private ProcessStages doNext;
		private TreeNode<T> next;
		private Iterator<TreeNode<T>> childrenCurNodeIter;
		private Iterator<TreeNode<T>> childrenSubNodeIter;

		public TreeNodeIter(TreeNode<T> treeNode) {
			this.treeNode = treeNode;
			this.doNext = ProcessStages.PROCESS_PARENT;
			this.childrenCurNodeIter = treeNode.children.iterator();
		}

		@Override
		public boolean hasNext() {

			if (this.doNext == ProcessStages.PROCESS_PARENT) {
				this.next = this.treeNode;
				this.doNext = ProcessStages.PROCESS_CHILD_CURRENT_NODE;
				return true;
			}

			if (this.doNext == ProcessStages.PROCESS_CHILD_CURRENT_NODE) {
				if (this.childrenCurNodeIter.hasNext()) {
					TreeNode<T> childDirect = this.childrenCurNodeIter.next();
					this.childrenSubNodeIter = childDirect.iterator();
					this.doNext = ProcessStages.PROCESS_CHILD_SUB_NODE;
					return hasNext();
				}

				else {
					this.doNext = null;
					return false;
				}
			}
			
			if (this.doNext == ProcessStages.PROCESS_CHILD_SUB_NODE) {
				if (this.childrenSubNodeIter.hasNext()) {
					this.next = this.childrenSubNodeIter.next();
					return true;
				}
				else {
					this.next = null;
					this.doNext = ProcessStages.PROCESS_CHILD_CURRENT_NODE;
					return hasNext();
				}
			}

			return false;
		}

		@Override
		public TreeNode<T> next() {
			return this.next;
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}

	}
	
}
