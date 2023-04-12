package net.sytes.botg.datatypes;

public class Parameter extends NamedValuePair {

	private Object targetValue;
	private Object upperLimit;
	private Object lowerLimit;
	private String unit;
	private String type;
	
	public Parameter() {		
	}
	
	public Parameter(String name, Object value) {
		super(name, value);
	}
	
	public Parameter(String name, Object value, Object lowerLimit, Object upperLimit, String unit, String type) {
		this(name, value);
		this.lowerLimit = lowerLimit;
		this.upperLimit = upperLimit;
		this.unit = unit;
		this.type = type;
	}
	
	public Object getTargetValue() {
		return this.targetValue;
	}

	public Object getUpperLimit() {
		return this.upperLimit;
	}

	public Object getLowerLimit() {
		return this.lowerLimit;
	}

	public String getUnit() {
		return this.unit;
	}

	public String getType() {
		return this.type;
	}
	
	@Override
	public String toString() {
		return this.getClass().getSimpleName() + "[" + this.name + "] = " + this.value;
	}
}