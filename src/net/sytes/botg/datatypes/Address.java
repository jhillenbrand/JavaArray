package net.sytes.botg.datatypes;

/**
 * class that stores data and information of address endpoint for a data node in data sources
 * @author botg@gmx.de
 */
public class Address {

	protected String source;
	protected String address;
	protected DataType dataType;
	protected String unit;
	protected String description;

	public Address() {
		this(new Builder());
	}
	
	private Address(Builder builder) {
		this.source = builder.source;
		this.address = builder.address;
		this.dataType = builder.dataType;
		this.unit = builder.unit;
		this.description = builder.description;
	}
	
	public static class Builder {
		
		private String source = "";
		private String address = "";
		private DataType dataType = DataType.DOUBLE;
		private String unit = "";
		private String description = "";

		public Builder source(String source) {
			this.source = source;
			return this;
		}

		public Builder address(String address) {
			this.address = address;
			return this;
		}

		public Builder dataType(DataType dataType) {
			this.dataType = dataType;
			return this;
		}
		
		public Builder dataType(String dataTypeStr) {
			this.dataType = DataType.convertToType(dataTypeStr);
			return this;
		}

		public Builder unit(String unit) {
			this.unit = unit;
			return this;
		}

		public Builder description(String description) {
			this.description = description;
			return this;
		}

		public Address build() {
			return new Address(this);
		}
	}
	
	/**
	 * makes a shallow copy of this object, as this.value can contain a object in {@code value}
	 * @return
	 */
	public Address clone() {
		Address newChannel = new Address.Builder()
				.source(this.source)
				.address(this.address)
				.dataType(this.dataType)
				.description(this.description)
				.unit(this.unit)
				.build();
		return newChannel;
	}
		
	public String toString() {
		StringBuilder sb = new StringBuilder();		
		sb.append(this.getClass().getSimpleName() + "\n");
		sb.append("\tsource = " + this.source + "\n");
		sb.append("\taddress = " + this.address + "\n");
		sb.append("dataType=" + this.dataType + "\n");
		sb.append("\tunit = " + this.unit + "\n");
		sb.append("\tdescription = " + this.description + "\n");
		return sb.toString();
	}

	public String getDescription() {
		return this.description;
	}
	
	public String getUnit() {
		return this.unit;
	}
	
	public String getAddress() {
		return this.address;
	}

	/**
	 * @return the dataType
	 */
	public DataType getDataType() {
		return this.dataType;
	}

	/**
	 * @return the source
	 */
	public String getSource() {
		return this.source;
	}
}
