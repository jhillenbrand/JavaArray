package net.sytes.botg.datatypes.buffers;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;

import net.sytes.botg.datastruct.ConfigOption;
import net.sytes.botg.datatypes.DataType;

/**
 * a synchronized FiFo buffer implementation, that allows to store Java {@code Object}s
 * <br>with a specified capacity
 * <br>the elements can be pulled and pushed element wise or as arrays 
 * @author jhillenbrand
 */
public class Buffer implements IBuffer {

	@ConfigOption
	protected String type = this.getClass().getName();
	
	@ConfigOption
	protected String id;
	
	@ConfigOption
	protected DataType dataType;
	
	@ConfigOption
	protected int capacity;
	
	@ConfigOption
	protected String description;
	
	@ConfigOption
	protected String unit;
	
	@ConfigOption
	protected Object[] initialValues;

	/** The buffered elements */
	protected Object[] elems;
	
	/** items index for next pull */
	protected int pullIndex;
	
	/** items index for push */
	protected int pushIndex;
	
	/** Number of elements in the queue */
	protected volatile int count;

	protected static final String VALUE_FIELD = "\"value\": ";
	protected static final String VALUES_FIELD = "\"values\": ";
	
	public static String uniqueId() {
		return Buffer.class.getSimpleName() + " [" + UUID.randomUUID() + "]";
	}
	
	public Buffer() {
		this(new Builder());
	}
	
	public Buffer(int capacity) {
		this(capacity, null);
	}
	
	public Buffer(int capacity, Object[] initialValues) {
		if (capacity <= 0) {
            throw new IllegalArgumentException("capacity must be greater than 0.");
    	}
        this.capacity = capacity;
		this.elems = new Object[capacity];
		this.initialValues = initialValues;
		if (this.initialValues != null) {
			
		}
	}
	
	private Buffer(Builder builder) {
		this(builder.capacity, builder.initialValues);
		this.id = builder.id;
		this.dataType = builder.dataType;
		this.capacity = builder.capacity;
		this.description = builder.description;
		this.unit = builder.unit;
	}

	public static class Builder {
		
		private String id = uniqueId();
		private DataType dataType = DataType.DOUBLE;
		private String description = null;
		private String unit = null;
		private int capacity = 1;
		private Object[] initialValues = null;

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
		
		public Builder initialValues(Object[] initialValues) {
			this.initialValues = initialValues;
			return this;
		}
		
		public Builder initialValues(List<Object> initialValues) {
			this.initialValues = initialValues.toArray();
			return this;
		}

		public Buffer build() {
			return new Buffer(this);
		}
	}

	@Override
	public synchronized void push(Object element) {
		 if (this.count == this.elems.length) {
         	this.dequeue();
         }
         this.enqueue(element);
	}
	
	public synchronized void push(Object[] elements) {
		this.enqueue(elements);
	}

	@Override
	public synchronized Object pull() {
		return this.dequeue();
	}

	@Override
	public synchronized Object newest() {
		int ni = this.pushIndex - 1;
		if (ni < 0) {
			ni = this.capacity() - 1;
		}
		Object element = this.elems[ni];
		return element;
	}

	@Override
	public synchronized Object oldest() {
		Object element = this.elems[this.pullIndex];    	
        return element;
	}

	@Override
	public synchronized Object[] toArray(boolean persistent) {
     	Object[] elements = this.copy();
        if (!persistent) { 
        	this.clear();
        }
        return elements;
	}
	
	@Override
	public synchronized Object[] toArray(boolean persistent, int numOfElements) {
		if (persistent) {
			return this.copy(numOfElements);
		} else {
			return this.dequeue(numOfElements);
		}
	}
		
	@Override
	public synchronized void clear() {
    	//this.elems = null; 
    	this.elems = new Object[this.capacity()];
    	this.pushIndex = 0;
    	this.pullIndex = 0;
        this.count = 0;
    }

	/**
	 * removes {@code numOfElements} Samples from {@code Buffer} in FiFo-order
	 * @param numOfElements
	 */
	public synchronized void clear(int numOfElements) {
		this.dequeue(numOfElements);
	}
	
	/**
	 * resets the {@code Buffer} to its original conditions
	 */
	public void reset() {
		this.clear();
		this.resize(this.capacity);
		if (this.initialValues != null) {
			 this.push(this.initialValues);
		}
	}

	/**
	 * resizes the current capacity of the {@code Buffer} to new {@code capacity}
	 */
	@Override
	public void resize(int capacity) {
		this.capacity = capacity;
		Object[] newElems = new Object[this.capacity];
		if (this.size() > 0) {
			System.arraycopy(this.elems, 0, newElems, 0, this.size());
		}
		this.elems = newElems;
	}
	
	@Override
	public int size() {
		return this.count;
	}
	
	@Override
	public int capacity() {
		if (this.capacity != this.elems.length) {
			this.resize(this.capacity);
		}
		return this.capacity;
	}
	
	@Override
	public int remainingCapacity() {
		return this.capacity() - this.size();
	}	
	
	/**
     * Inserts element at tail of buffer
     */
    protected void enqueue(Object element) {
        this.elems[this.pushIndex] = element;
    	if (++this.pushIndex == this.elems.length) {
        	this.pushIndex = 0;
        }
        this.count++;
    }
    
    /**
     * Inserts {@code elements} as new elements into this {@code Buffer}
     * <br>must be called in a synchronized method
     * @param elements
     */
    protected void enqueue(Object[] elements) {
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
			this.pushIndex = this.pushIndex + m;
			if (this.pushIndex == c) {
	        	this.pushIndex = 0;	        	
	        }
			this.count = m;
		
		// b) remaining capacity is enough to store new elements 
		} else if (m <= rm) {
			System.arraycopy(elements, 0, this.elems, this.pushIndex, m);
			this.pushIndex = this.pushIndex + m;
			if (this.pushIndex == c) {
	        	this.pushIndex = 0;
	        }
			this.count = this.count + m;
		
		// c) remaining capacity is not enough to store new elements, but the capacity is enough
		} else if (m >= rm && m <= c) {
			if (c - this.pushIndex >= m) {
				System.arraycopy(elements, 0, this.elems, this.pushIndex, m);
				this.pushIndex = this.pushIndex + m;
				this.pullIndex = this.pushIndex;	 					
			} else {
				int rem1 = c - this.pushIndex;
				int rem2 = m - rem1;
				System.arraycopy(elements, 0, this.elems, this.pushIndex, rem1);
				System.arraycopy(elements, rem1, this.elems, 0, rem2);
				this.pullIndex = rem2;
				this.pushIndex = rem2;
			}
			this.count = c;
			
		// d) new elements are too much for capacity
		} else if (m > c) {
			int s = m - c;
			System.arraycopy(elements, s, this.elems, 0, c);
			this.pushIndex = 0;
			this.pullIndex = 0;
			this.count = c;
			
		// e) case should not happen	
		} else {
			throw new IllegalStateException("Unknown State for adding multiple Buffer elements  ");
		}		
		
    }

    /**
     * removes head of buffer
     */
    protected Object dequeue() {
    	Object element = this.elems[this.pullIndex];    	
    	this.elems[this.pullIndex] = null;
    	if (++this.pullIndex == this.elems.length) {
            this.pullIndex = 0;
        }
        this.count--;
        return element;
    }
    
    /**
     * removes {@code numOfElements} elements from the {@code Buffer} into an array
     * <br>must be called within a synchronized method
     * @param numOfElements
     * @return
     */
    protected Object[] dequeue(int numOfElements) {
    	int count = this.count;
		if (numOfElements > count) {
			numOfElements = count;
		}
		Object[] elements = this.copy(numOfElements);
		Object[] delElements = new Object[numOfElements];
		
		int n = this.elems.length - this.pullIndex;
		if (numOfElements <= n) {
            System.arraycopy(delElements, 0, this.elems, this.pullIndex, numOfElements);
        } else {
            System.arraycopy(delElements, 0, this.elems, this.elems.length - n, n);
            System.arraycopy(delElements, n, this.elems, 0, numOfElements - n);
        }
		// correct pullIndex and count
        this.pullIndex = this.pullIndex + numOfElements;
        int d = this.pullIndex - this.elems.length;
        if (d == 0) {
            this.pullIndex = 0;
        } else if (d > 0){
        	this.pullIndex = d;
        } 
        this.count = count - numOfElements;
        if (this.count == 0) {
        	this.pullIndex = 0;
        	this.pushIndex = 0;
        }
		return elements;    	
    }
    
    /**
     * copies all elements from {@code Buffer} into array
     * @return
     */
    protected Object[] copy() {
    	int count = this.count;
		Object[] elements = new Object[count];
        int n = this.elems.length - this.pullIndex;
        if (count <= n) {
            System.arraycopy(this.elems, this.pullIndex, elements, 0, count);
        } else {
            System.arraycopy(this.elems, this.pullIndex, elements, 0, n);
            System.arraycopy(this.elems, 0, elements, n, count - n);
        }
        return elements;
    }
    
    /**
     * copies {@code numOfElements} from {@code Buffer} to array
     * <br>if {@code numOfElements} exceeds the number of available elements, only the available elements are copied
     * <br>must be called in a synchronized method
     * @param numOfElements
     * @return
     */
    protected Object[] copy(int numOfElements) {
    	int count = this.count;
		if (numOfElements > count) {
			numOfElements = count;
		}
    	Object[] elements = new Object[numOfElements];
		int n = this.elems.length - this.pullIndex;
		if (numOfElements <= n) {
            System.arraycopy(this.elems, this.pullIndex, elements, 0, numOfElements);
        } else {
            System.arraycopy(this.elems, this.pullIndex, elements, 0, n);
            System.arraycopy(this.elems, 0, elements, n, numOfElements - n);
        }
		return elements;
    }
	
    /**
     * method blocks until buffer has reached specified {@code size},
     * <br>the condition is check every {@code updateRate} ms
     * @param size	
     * @param updateRate [ms]
     */
    public void waitForSize(int size, int updateRate) {
    	if (size > this.capacity()) {
    		throw new IllegalArgumentException("size must be smaller than capacity of " + this.getClass().getSimpleName() + " (" + this.capacity() + ")");
    	}
    	while(this.size() < size) {
    		try {
				Thread.sleep(updateRate);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
    	}
    }
	
    /**
     * returns this {@code Buffer}s config as a {@code Map<String, Object>}
      * @return
     */
	public Map<String, Object> toMap(){
		Map<String, Object> map = new LinkedHashMap<String, Object>();
		map.put("type", this.type);
		map.put("id", this.id);
		map.put("capacity", this.capacity());
		map.put("size", this.size());
		map.put("dataType", this.dataType.toString());
		map.put("unit", this.unit);
		map.put("description", this.description);
		if (this.initialValues != null) {
			map.put("initialValues", Arrays.asList(this.initialValues));
		}
		return map;
	}
    
	/**
     * returns a String describing content and config
     */
	@Override
    public String toString() {
    	return this.toMap().toString();    	
    }
	
	/**
     * returns a String describing content, config and data (if {@code withData} is set to true)
     */
    public String toString(boolean withData) {
    	if (withData) {
	    	StringBuilder sb = new StringBuilder();
	    	sb.append(this.getClass().getSimpleName()).append(" id=[").append(this.id).append("]");    	
	    	sb.append("\n\tElements:" + this.size() + "/" + this.capacity());
	    	sb.append("\n\tDataType: " + this.dataType);
	    	sb.append("\n\t" + Arrays.toString(unsynchronizedToArray()));
			return sb.toString();
    	} else {
    		return this.toString();
    	}    	
    }
    
    /**
     * views the current buffer elements without synchronization [TODO MIGHT BE UNSTABLE]
     * @return
     */
    protected Object[] unsynchronizedToArray() {
		Object[] elements;
     	final int count = this.count;
     	final int pullIndex = this.pullIndex;
     	elements = new Object[count];
        int n = this.elems.length - pullIndex;
        if (count <= n) {
            System.arraycopy(this.elems, pullIndex, elements, 0, count);
        	} else {
            System.arraycopy(this.elems, pullIndex, elements, 0, n);
            System.arraycopy(this.elems, 0, elements, n, count - n);
        }
        return elements;
	}

	@Override
	public String jsonValue() {
		StringBuilder sb = new StringBuilder();
		sb.append("{\n\t").append(VALUE_FIELD).append(this.newest().toString()).append("\n}");		
		return sb.toString();
	}

	@Override
	public String jsonValues() {
		StringBuilder sb = new StringBuilder();
		sb.append("{\n\t").append(VALUES_FIELD).append(Arrays.toString(this.toArray(true))).append("\n}");		
		return sb.toString();
	}

	public String getType() {
		return this.type;
	}

	public String getId() {
		return this.id;
	}

	public DataType getDataType() {
		return this.dataType;
	}

	public String getDescription() {
		return this.description;
	}

	public String getUnit() {
		return this.unit;
	}
}
