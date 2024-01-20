package net.sytes.botg.datastruct;

import java.util.UUID;

public class NamedValuePair implements ValueInterface {
	
	protected String name;
	
	/**
	 * property that stores the data of this class object
	 * http://tutorials.jenkov.com/java-concurrency/volatile.html
	 * The Java volatile keyword guarantees visibility of changes to variables across threads.
	 */
	protected volatile Object value;
		
	public NamedValuePair() {
		this(NamedValuePair.class.getSimpleName() + "[" + UUID.randomUUID().toString() + "]", null);
	}
	
	public NamedValuePair(String name, Object value) {
		this.name = name;
		this.value = value;
	}
	
	public String getName() {
		return this.name;
	}

	@Override
	public Object getValue() {
		return this.value;
	}

	@Override
	public void setValue(Object value) {
		this.value = value;
	}

	@Override
	public void setName(String name) {
		this.name = name;
	}
}
