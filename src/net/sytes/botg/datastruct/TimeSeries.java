package net.sytes.botg.datastruct;

import java.util.Arrays;

import net.sytes.botg.datatypes.buffers.TimedBuffer;

public class TimeSeries {
	
	private static String VALUES_IDENTIFIER = "values";
	private static String TIMESTAMPS_IDENTIFIER = "timestamps";
	
	private long[] timestamps = null;
	private Object[] values = null;
	
	/**
	 * create an empty {@code TimeSeries} of size {@code n}
	 * <br>this object can later be populated with the {@code setSample(int i, long time, Object value)} Method
	 * @param n
	 */
	public TimeSeries(int n) {
		this.timestamps = new long[n];
		this.values = new Object[n];
	}
	
	public TimeSeries(long[] timestamps, Object[] values) {
		this.timestamps = timestamps;
		this.values = values;
	}

	public TimeSeries(Long[] t, Object[] v) {
		this(toPrimitive(t), v);
	}
	
	public <E> TimeSeries(Sample<E>[] samples) {
		this.timestamps = new long[samples.length];
		this.values = new Object[samples.length];
		int s = 0;
		for (Sample<E> sample : samples) {
			this.timestamps[s] = sample.time();
			this.values[s] = sample.value();
			++s;
		}
	}
	
	/**
	 * sets the sample of the {@code TimeSeries} at index {@code i}
	 * @param i
	 * @param value
	 * @param time
	 */
	public void setSample(int i, Object value, long time) {
		this.values[i] = value;
		this.timestamps[i] = time;
	}
	
	/**
	 * returns the number of samples in this {@code TimeSeries}
	 * @return
	 */
	public int size() {
		return this.values.length;
	}

	public long[] getTimestamps() {
		return this.timestamps;
	}

	public Object[] getValues() {
		return this.values;
	}
	
	/**
	 * 
	 * @param timedBuffer
	 * @return
	 */
	public static TimeSeries fromBuffer(TimedBuffer timedBuffer) {
		return timedBuffer.toTimeSeries(true);
	}
	 
	/**
	 * inserts the{@code TimeSeries} data into the {@code timedBuffer}
	 * @param timedBuffer
	 */
	public void toBuffer(TimedBuffer timedBuffer) {
		timedBuffer.push(this.timestamps, this.values);
	}
	
	private static long[] toPrimitive(final Long[] array) {
        if (array == null) {
            return null;
        } else if (array.length == 0) {
            return new long[]{};
        }
        final long[] result = new long[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].longValue();
        }
        return result;
    }
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(this.getClass().getSimpleName())
			.append(":\n\t")
			.append(TIMESTAMPS_IDENTIFIER)
			.append(" = ")
			.append(Arrays.toString(this.timestamps))
			.append("\n\t ")
			.append(VALUES_IDENTIFIER)
			.append(" = ")
			.append(Arrays.toString(this.values))
			.append("\n");
		
		return sb.toString();
	}
	
}
