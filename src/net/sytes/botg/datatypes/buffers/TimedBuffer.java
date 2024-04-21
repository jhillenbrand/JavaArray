package net.sytes.botg.datatypes.buffers;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.UUID;

import net.sytes.botg.datastruct.Sample;
import net.sytes.botg.datastruct.TimeSeries;
import net.sytes.botg.datatypes.DataType;

public class TimedBuffer extends Buffer implements ITimedBuffer {
	
	protected long[] times;
	
	protected static final String TIME_FIELD = "\"time\": ";
	protected static final String TIMES_FIELD = "\"times\": ";
		
    public static Map<String, TimedBuffer> toMap(TimedBuffer buffer){
    	Map<String, TimedBuffer> buffers = new LinkedHashMap<String, TimedBuffer>();
    	buffers.put(buffer.getId(), buffer);
    	return buffers;
    }
	
	public TimedBuffer() {
		this(new Builder());
	}
	
	public TimedBuffer(int capacity) {
		if (capacity <= 0) {
            throw new IllegalArgumentException("Capacity must be greater than 0.");
    	}
		this.capacity = capacity;
		this.elems = new Object[capacity];
		this.times = new long[capacity];
	}
	
	private TimedBuffer(Builder builder) {
		this(builder.capacity);
		this.id = builder.id;
		this.dataType = builder.dataType;
		this.capacity = builder.capacity;
		this.description = builder.description;
		this.unit = builder.unit;
	}

	public static class Builder {
		
		private String id = TimedBuffer.class.getSimpleName() + " [" + UUID.randomUUID().toString() + "]";
		private DataType dataType = DataType.DOUBLE;
		private String description = null;
		private String unit = null;
		private int capacity = 1;

		public Builder id(String id) {
			this.id = id;
			return this;
		}

		public Builder dataType(DataType dataType) {
			this.dataType = dataType;
			return this;
		}
		
		public Builder description(String description) {
			this.description = description;
			return this;
		}

		public Builder unit(String unit) {
			this.unit = unit;
			return this;
		}
		
		public Builder capacity(int capacity) {
			this.capacity = capacity;
			return this;
		}

		public TimedBuffer build() {
			return new TimedBuffer(this);
		}
	}

	@Override
	public synchronized Sample<Object> pull() {
		return this.dequeue();
	}
	
	@Override
	public synchronized void push(long time, Object element) {
		 if (this.count == this.elems.length) {
         	this.dequeue();
         }
         this.enqueue(time, element);
	}
	
	@Override
	public synchronized void push(long[] times, Object[] elements) {
		this.enqueue(times, elements);
	}

	@Override
	public synchronized Sample<Object> newest() {
		int ni = this.pushIndex - 1;
		if (ni < 0) {
			ni = this.capacity() - 1;
		}
		Object element = this.elems[ni];
		long time = this.times[ni];
		if (element != null) {
			return new Sample<Object>(time, element);
		} else {
			return null;
		}
	}

	@Override
	public synchronized Sample<Object> oldest() {
		Object element = this.elems[this.pullIndex]; 
		long time = this.times[this.pullIndex];
		if (element != null) {
			return new Sample<Object>(time, element);
		} else {
			return null;
		}
	}

	@Override
	public synchronized TimeSeries toTimeSeries(boolean persistent) {
		TimeSeries ts = this.copyTimeSeries();
        if (!persistent) { 
        	this.clear();
        }
        return ts;
	}

	@Override
	public synchronized TimeSeries toTimeSeries(boolean persistent, int n) {
		return this.toTimeSeries(persistent, n, false);
	}
	
	@Override
	public synchronized TimeSeries toTimeSeries(boolean persistent, int n, boolean blocking) {
		if (!blocking) {
			if (persistent) {
				return this.copyTimeSeries(n);
			} else {
				return this.dequeueTimeSeries(n);
			}
		} else {
			this.waitForSize(n, 0);
			if (persistent) {
				return this.copyTimeSeries(n);
			} else {
				return this.dequeueTimeSeries(n);
			}
		}
	}
	
	@Override
	public synchronized void clear() {
    	super.clear();
		this.times = new long[this.capacity()];
    }
	
	/**
	 * removes {@code numOfSamples} Samples from {@code TimedBuffer} in FiFo-order
	 * @param numOfSamples
	 */
	public synchronized void clear(int numOfSamples) {
		this.dequeueTimeSeries(numOfSamples);
	}
	
	/**
     * Inserts {@code time} and {@code element} at tail of {@code TimedBuffer}
     */
    protected void enqueue(long time, Object element) {
        this.elems[this.pushIndex] = element;
        this.times[this.pushIndex] = time; 
    	if (++this.pushIndex == this.elems.length) {
        	this.pushIndex = 0;
        }
        this.count++;
    }
      
    /**
     * Inserts {@code times} and {@code elements} as new elements into this {@code TimedBuffer}
     * <br>must be called within synchronized method
     * @param times
     * @param elements
     */
    protected void enqueue(long[] times, Object[] elements) {
    	if (times == null && elements != null) {
    		super.enqueue(elements);
    		return;
    	}
    	if (times.length != elements.length) {
    		StringBuilder sb = new StringBuilder();
    		sb.append("Number of elements and times must be the same (")
    			.append(times.length)
    			.append(" <> ")
    			.append(elements.length)
    			.append(")");
    		throw new IllegalArgumentException(sb.toString());
    	}
    	
    	// differentiate between cases
		
    	final int c = this.capacity();
    	final int rm = this.remainingCapacity();
		final int m = elements.length;
		
		// check if pushindex has reached end
		//if (this.pushIndex == c) {
        //	this.pushIndex = 0;
        //}
		
		// a) new elements are smaller or equal to buffer capacity and buffer is empty
		if (m <= c && this.count == 0) {
			System.arraycopy(elements, 0, this.elems, 0, m);
			System.arraycopy(times, 0, this.times, 0, m);
			this.pushIndex = this.pushIndex + m;
			if (this.pushIndex == c) {
	        	this.pushIndex = 0;
	        }
			this.count = m;
		
		// b) remaining capacity is enough to store new elements 
		} else if (m <= rm) {
			System.arraycopy(elements, 0, this.elems, this.pushIndex, m);
			System.arraycopy(times, 0, this.times, this.pushIndex, m);
			this.pushIndex = this.pushIndex + m;
			if (this.pushIndex == c) {
	        	this.pushIndex = 0;
	        }
			this.count = this.count + m;
					
		// c) remaining capacity is not enough to store new elements, but the capacity is enough
		} else if (m >= rm && m <= c) {
			if (c - this.pushIndex >= m) {
				System.arraycopy(elements, 0, this.elems, this.pushIndex, m);
				System.arraycopy(times, 0, this.times, this.pushIndex, m);
				this.pushIndex = this.pushIndex + m;
				this.pullIndex = this.pushIndex;	 					
			} else {
				int rem1 = c - this.pushIndex;
				int rem2 = m - rem1;
				System.arraycopy(elements, 0, this.elems, this.pushIndex, rem1);
				System.arraycopy(elements, rem1, this.elems, 0, rem2);
				System.arraycopy(times, 0, this.times, this.pushIndex, rem1);
				System.arraycopy(times, rem1, this.times, 0, rem2);
				this.pullIndex = rem2;
				this.pushIndex = rem2;
			}
			this.count = c;
			
		// d) new elements are too much for capacity
		} else if (m > c) {
			int s = m - c;
			System.arraycopy(elements, s, this.elems, 0, c);
			System.arraycopy(times, s, this.times, 0, c);
			this.pushIndex = 0;
			this.pullIndex = 0;
			this.count = c;
			
		// e) case should not happen	
		} else {
			throw new IllegalStateException("Unknown State for adding multiple elements into " + this.getClass().getSimpleName());
		}	
    }

    /**
     * removes head of buffer
     */
    protected Sample<Object> dequeue() {
    	long time = this.times[this.pullIndex];
    	Object element = this.elems[this.pullIndex];
    	if (element != null) {
	    	this.elems[this.pullIndex] = null;
	    	this.times[this.pullIndex] = 0; 
	    	if (++this.pullIndex == this.elems.length) {
	            this.pullIndex = 0;
	        }
	        this.count--;
	        return new Sample<Object>(time, element);
    	} else {
    		return null;
    	}
    }
    
    /**
     * removes {@code numOfSamples} samples from the {@code TimedBuffer} into an {@code TimeSeries}
     * <br>must be called within a synchronized method
     * @param numOfSamples
     * @return
     */
    protected TimeSeries dequeueTimeSeries(int numOfSamples) {
    	int count = this.count;
    	if (count == 0) {
    		return null;
    	}
		if (numOfSamples > count) {
			numOfSamples = count;
		}
		TimeSeries ts = this.copyTimeSeries(numOfSamples);
		Object[] delElements = new Object[numOfSamples];
		long[] delTimes = new long[numOfSamples];
		
		int n = this.elems.length - this.pullIndex;
		if (numOfSamples <= n) {
            System.arraycopy(delElements, 0, this.elems, this.pullIndex, numOfSamples);
            System.arraycopy(delTimes, 0, this.times, this.pullIndex, numOfSamples);
        } else {
            System.arraycopy(delElements, 0, this.elems, this.elems.length - n, n);
            System.arraycopy(delElements, n - 1, this.elems, 0, numOfSamples - n);
            System.arraycopy(delTimes, 0, this.times, this.times.length - n, n);
            System.arraycopy(delTimes, n - 1, this.times, 0, numOfSamples - n);
        }
		// correct pullIndex and count
        this.pullIndex = this.pullIndex + numOfSamples;
        int d = this.pullIndex - this.elems.length;
        if (d == 0) {
            this.pullIndex = 0;
        } else if (d > 0){
        	this.pullIndex = d;
        } 
        this.count = count - numOfSamples;
        if (this.count == 0) {
        	this.pullIndex = 0;
        	this.pushIndex = 0;
        }
    	return ts;    	
    }
    
    /**
     * copies this {@code TimedBuffer} content as {@code TimeSeries}
     * <br>must be called in a synchronized method
     * @return
     */
    protected TimeSeries copyTimeSeries() {
    	long[] times;
		Object[] elements;
     	final int count = this.count;
     	if (count == 0) {
     		return null;
     	}
     	elements = new Object[count];
     	times = new  long[count];
        int n = this.elems.length - this.pullIndex;
        if (count <= n) {
            System.arraycopy(this.elems, this.pullIndex, elements, 0, count);
            System.arraycopy(this.times, this.pullIndex, times, 0, count);
    	} else {
            System.arraycopy(this.elems, this.pullIndex, elements, 0, n);
            System.arraycopy(this.elems, 0, elements, n, count - n);
            System.arraycopy(this.times, this.pullIndex, times, 0, n);
            System.arraycopy(this.times, 0, times, n, count - n);
        }
        return new TimeSeries(times, elements);
    }
    
    /**
     * copies {@code numOfSamples} from {@code TimedBuffer} to {@code TimeSeries}
     * <br>if {@code numOfSamples} exceeds the number of available samples, only the available samples are copied
     * <br>must be called in a synchronized method
     * @param samples
     * @return
     */
    protected TimeSeries copyTimeSeries(int numOfSamples) {
    	int count = this.count;
    	if (count == 0) {
    		return null;
    	}
		if (numOfSamples > count) {
			numOfSamples = count;
		}
		long[] times = new long[numOfSamples];
    	Object[] elements = new Object[numOfSamples];
		int n = this.elems.length - this.pullIndex;
		if (numOfSamples <= n) {
			System.arraycopy(this.times, this.pullIndex, times, 0, numOfSamples);
            System.arraycopy(this.elems, this.pullIndex, elements, 0, numOfSamples);            
        } else {
        	System.arraycopy(this.times, this.pullIndex, times, 0, n);
            System.arraycopy(this.times, 0, times, n, numOfSamples - n);
            System.arraycopy(this.elems, this.pullIndex, elements, 0, n);
            System.arraycopy(this.elems, 0, elements, n, numOfSamples - n);
        }
		/*
		if (elements[0] == null) {
			System.out.println("should not happen");
		}
		*/
		return new TimeSeries(times, elements);
    }
	
	/**
    * returns a String describing content and config
    */
	@Override
	public String toString() {
		return super.toString();     	
	}
	
	/**
    * returns a String describing content, config and data (if {@code withData} is set to true)
    */
	@Override
	public String toString(boolean withData) {
		if (withData) {
			StringBuilder sb = new StringBuilder();
			sb.append(this.getClass().getSimpleName()).append(" id=[").append(this.id).append("]");    	
			sb.append("\n\tElements:" + this.size() + "/" + this.capacity());
			sb.append("\n\tDataType: " + this.dataType);
			sb.append("\n\t" + unsynchronizedToTimeSeries());
			return sb.toString();
    	} else {
    		return this.toString();
    	}    	
	}
	
	/**
	 *  returns the newest buffer value as json string
	 *  <br>
	 *  <br>Example:
	 *  <br>{
	 *  <br>&nbsp;"time": 1603403432210,
	 *  <br>&nbsp;"value": 10.0 
	 *  <br>}
	 * @return
	 */
	@Override
	public String jsonValue() {
		StringBuilder sb = new StringBuilder();
		Sample<Object> sample = this.newest();
		sb.append("{\n\t").append(TIME_FIELD).append(sample.time());
		sb.append("\n\t").append(VALUE_FIELD).append(sample.value().toString());
		sb.append("\n}");		
		return sb.toString();
	}

	/**
	 *  returns the newest buffer value as json string
	 *  <br>
	 *  <br>Example:
	 *  <br>{
	 *  <br>&nbsp;"times": [1603403432210, 1603403432310, 1603403432314],
	 *  <br>&nbsp;"values": [10.0, 9.23, 8.12] 
	 *  <br>}
	 * @return
	 */
	@Override
	public String jsonValues() {
		StringBuilder sb = new StringBuilder();
		TimeSeries ts = this.toTimeSeries(true);
		sb.append("{\n\t").append(TIMES_FIELD).append(Arrays.toString(ts.getTimestamps()).replace("{", "[").replace("}", "]"));
		sb.append("\n\t").append(VALUES_FIELD).append(Arrays.toString(ts.getValues()).replace("{", "[").replace("}", "]"));
		sb.append("\n}");		
		return sb.toString();
	}
   
	protected TimeSeries unsynchronizedToTimeSeries() {
		long[] times;
		Object[] elements;
     	final int count = this.count;
     	final int pullIndex = this.pullIndex;
     	elements = new Object[count];
     	times = new  long[count];
        int n = this.elems.length - pullIndex;        
        if (count <= n) {
            System.arraycopy(this.elems, pullIndex, elements, 0, count);
            System.arraycopy(this.times, pullIndex, times, 0, count);
    	} else { 
            System.arraycopy(this.elems, pullIndex, elements, 0, n);
            System.arraycopy(this.elems, 0, elements, n, count - n);
            System.arraycopy(this.times, pullIndex, times, 0, n);
            System.arraycopy(this.times, 0, times, n, count - n);
        }
        return new TimeSeries(times, elements);
	}
}
