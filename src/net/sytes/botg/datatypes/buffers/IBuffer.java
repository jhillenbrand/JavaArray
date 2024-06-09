package net.sytes.botg.datatypes.buffers;

public interface IBuffer {

	/**
	 * inserts a new element to buffer, if buffer was full, it discards the oldest before adding new element
	 * @param element
	 */
	public void push(Object element);
	
	/**
	 * inserts new elements to buffer, if buffer becomes full during inserting, old elements are removed accordingly, when new ones are added
	 * @param elements
	 */
	public void push(Object[] elements);
	
	/**
	 * returns and removes the oldest (head) element of the buffer
	 * @return
	 */
	public Object pull();
	
	/**
	 * returns the newest (tail) element of the buffer
	 * @return
	 */
	public Object newest();
	
	/**
	 * returns the oldest (head) element of the buffer
	 * @return
	 */
	public Object oldest();
	
	/**
	 * returns all elements in the buffer in FiFo order, if {@code persistent} is {@code true}, then the buffer is {@code  clear}ed afterwards
	 * @param persistent
	 * @return
	 */
	public Object[] toArray(boolean persistent);
	
	/**
	 * returns {@code numOfElements} elements of {@code Buffer} 
	 * @param persistent
	 * @param numOfElements
	 * @return
	 */
	public Object[] toArray(boolean persistent, int numOfElements);
	
	/**
	 * removes all elements from the buffer
	 */
	public void clear();
	
	/**
	 * resizes the {@code capacity} of this {@code Buffer}	
	 * @param capacity
	 */
	public void resize(int capacity);
	
	/**
	 * returns the current size (number of elements in the buffer)
	 * @return
	 */
	public int size();
	
	/**
	 * returns the maximum size of the buffer
	 * @return
	 */
	public int capacity();
	
	/**
	 * returns the remaining number of elements, that can be added
	 * @return
	 */
	public int remainingCapacity();
	
	/**
	 *  returns the newest buffer value as json string
	 *  <br>
	 *  <br>Example:
	 *  <br>{
	 *  <br>&nbsp;"value": 10.0
	 *  <br>}
	 * @return
	 */
	public String jsonValue();
	
	/**
	 * returns all buffer values as json string
	 * <br>
	 * <br>Example:
	 * <br>{
	 * <br>&nbsp;"values": [10.0, 9.1, 87.1, 0.10]
	 * <br>}
	 * @return
	 */
	public String jsonValues();
	
}
