package net.sytes.botg.datatypes;

public class ConfigPair {

	private String attribute = null;
	private String value = null;
	
	public ConfigPair() {
	}
	
	public ConfigPair(String attribute, String value) {
		this.attribute = attribute;
		this.value = value;
	}
	
	public String getValue() {
		return this.value;
	}
	
	public void setValue(Object value) {
		this.value = (String) value;
	}

	public String getAttribute() {
		return this.attribute;
	}

	public void setAttribute(String name) {
		this.attribute = name;
	}
	
	@Override
	public String toString() {
		return this.attribute + " = " + this.value; 
	}

}
