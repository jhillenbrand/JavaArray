package net.sytes.botg.datatypes.buffers;

import net.sytes.botg.datastruct.Sample;
import net.sytes.botg.datastruct.TimeSeries;

public interface ITimedBuffer {

	/**
	 * returns the oldest (head) Sample {time, value} from the buffer
	 * @return
	 */
	public Sample<Object> pull();
	
	/**
	 * inserts a new element and time to the buffer (tail), if buffer is full, then it discards the oldest element before adding new element
	 * @param t
	 * @param element
	 */
	public void push(long time, Object element);
	
	/**
	 * inserts elements and times to the buffer (tail), if buffer is full, then it discards enough elements to add the new elements
	 * @param time
	 * @param element
	 */
	public void push(long[] times, Object[] elements);
	
	/**
	 * returns the latest added {@code Sample} from {@code TimedBuffer}
	 * @return
	 */
	public Sample<Object> newest();
	
	/**
	 * returns the earliest added {@code Sample} from {@code TimedBuffer}
	 * @return
	 */
	public Sample<Object> oldest();
	
	/**
	 * returns the entire {@code TimedBuffer} as {@code TimeSeries}
	 * <br>if persistent is set to false, the {@code Sample}s in the {@code TimedBuffer} are removed
	 * @param persistent
	 * @return
	 */
	public TimeSeries toTimeSeries(boolean persistent);
	
	/**
	 * returns {@code numOfSamples} from {@code TimedBuffer} as {@code TimeSeries}
	 * <br>if persistent is set to false, the extracted {@code Sample}s in the {@code TimedBuffer} are removed
	 * @param persistent
	 * @param numOfSamples
	 * @return
	 */
	public TimeSeries toTimeSeries(boolean persistent, int numOfSamples);

	/**
	 * returns the content of this {@code TimedBuffer} as {@code TimeSeries}
	 * <br>if persistent is set to false, the elements of the buffer requested are being removed
	 * <br>only {@code samples} of the buffer are retrieved
	 * <br>if {@code blocking} is set to true, this methods blocks until enough elements are in {@code TimedBuffer} to retrieve {@code n} samples
	 * @param persistent
	 * @param n
	 * @param blocking
	 * @return
	 */
	TimeSeries toTimeSeries(boolean persistent, int n, boolean blocking);
		
}
