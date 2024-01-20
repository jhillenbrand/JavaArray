package net.sytes.botg.datastruct;

/**
 * representing a value pair of {time, value}
 * @author hillenbrand *
 * @param <T>
 */
public class Sample<T> {
	
	private final long time;
	private final T value;
	
	public Sample(T value) {
		this.time = 0;
		this.value = value;
	}
	
    public Sample(long time, T value) {
    	this.time = time;
    	this.value = value;
    }
    
    public long time() {
		return this.time;
	}

	public T value() {
		return this.value;
	}
    
    @Override
    public String toString() {
    	return "(" + this.time + ", " + this.value.toString() + ")"; 
    }
}
