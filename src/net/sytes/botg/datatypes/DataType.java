package net.sytes.botg.datatypes;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public enum DataType {
	
	DOUBLE(Double.BYTES), INT(Integer.BYTES), FLOAT(Float.BYTES), LONG(Long.BYTES), CHAR(Character.BYTES), STRING, BYTE, BOOLEAN(1), SHORT(Short.BYTES), OBJECT;

	private int bytes = -1;

	private DataType(int bytes) {
		this.bytes = bytes;
	}
	
	private DataType(){
		this(-1);
	}
	
	/**
	 * The provided code is a Java method named stringToDataType that takes a String argument named dataTypeStr and returns an instance of the DataType enum. The method is marked as public and static, which means it can be called from other classes without requiring an instance of the class.
	 * The method first checks if the dataTypeStr argument is null, and if so, it returns null indicating that the provided argument is not valid.
	 * The method then converts the provided dataTypeStr argument to lowercase using the toLowerCase method of the String class to ensure that case sensitivity is not an issue in the comparison.
	 * Next, the method uses a series of if-else statements to determine which DataType enum instance should be returned based on the dataTypeStr argument. The if-else statements check if the str variable matches any of the known data types and returns the corresponding enum value.
	 * Finally, if the provided dataTypeStr argument does not match any of the known data types, the method returns null.
	 * Overall, the provided code appears to be a useful utility method for converting a string representation of a data type into its corresponding DataType enum value.
	 * @param dataTypeStr
	 * @return
	 */
	public static DataType stringToDataType(String dataTypeStr) {
		if (dataTypeStr == null) {
			return null;
		}
		String str = dataTypeStr.toLowerCase();
		if (str.contentEquals("int64") || str.contentEquals("long") || str.contentEquals("lint")) {
			return LONG;
		} else if (str.contentEquals("int32") || str.contentEquals("int") || str.contentEquals("integer") || str.contentEquals("dint")) {
			return INT;
		} else if (str.contentEquals("int16") || str.contentEquals("short") || str.contentEquals("shortint") || str.contentEquals("smallint")) {
			return SHORT;
		} else if (str.contentEquals("double") || str.contentEquals("lreal")) {
			return DOUBLE;
		} else if (str.contentEquals("string") || str.contentEquals("java.lang.string")) {
			return STRING;
		} else if (str.contentEquals("char") || str.contentEquals("char[]")) {
			return CHAR;
		} else if (str.contentEquals("bool") || str.contentEquals("boolean") || str.contentEquals("bit")) {
			return BOOLEAN;
		} else if (str.contentEquals("float") || str.contentEquals("real")) {
			return FLOAT;
		} else if (str.contentEquals("byte") || str.contentEquals("bytes") || str.contentEquals("byte[]") || str.contentEquals("int8")) {
			return BYTE;
		} else if (str.contentEquals("object") || str.contentEquals("obj") || str.contentEquals("object[]")) {
			return OBJECT;
		} else {
			return null;
		}
	}
	
	/**
	 * converts the passed object to {@code targetType} (by determining the class of  {@code obj} beforehand)
	 * @param obj
	 * @param targetType
	 * @return
	 */
	public static Object cast(Object obj, DataType targetType) {
		return cast(obj, of(obj.getClass()), targetType);
	}
	
	/**
	 * converts the passed object from  {@code sourceType} to {@code targetType}
	 * @param obj
	 * @param sourceType
	 * @param targetType
	 * @return
	 */
	public static Object cast(Object obj, DataType sourceType, DataType targetType) throws NumberFormatException {
		String s = null;
		switch(sourceType) {
			case BOOLEAN:
				switch (targetType) {
					case BOOLEAN:
						return (boolean) obj;
						
					case BYTE:
						return (byte) (((boolean) obj) ? 1 : 0);
						
					case CHAR:
						return (char) (((boolean) obj) ? '1' : '0');
						
					case DOUBLE:
						return (double) (((boolean) obj) ? 1.0 : 0.0);
						
					case FLOAT:
						return (float) (((boolean) obj) ? 1.0f : 0.0f);
						
					case INT:
						return (int) (((boolean) obj) ? 1 : 0);
						
					case LONG:
						return (long) (((boolean) obj) ? 1L : 0L);
						 
					case SHORT:
						short sh;
						sh = (short) (((boolean) obj) ? 1 : 0);
						return sh;
						
					case STRING:
						return ((boolean) obj) ? "true" : "false";
						
					case OBJECT:
						return obj;
						
					default:
						return null;				
				}
				
			case BYTE:
				switch (targetType) {
					case BOOLEAN:
						return (boolean) (((byte) obj) == 1 ? true : false);
						
					case BYTE:
						return (byte) obj;
						
					case CHAR:
						return (char) ((byte) obj);
						
					case DOUBLE:
						return (double) ((byte) obj);
						
					case FLOAT:
						return (float) ((byte) obj);
						
					case INT:
						return (int) ((byte) obj);
						
					case LONG:
						return (long) ((byte) obj);					
						
					case SHORT:
						return (short) ((byte) obj);
												
					case STRING:
						return new String(new byte[]{(byte) obj});
					
					case OBJECT:
						return obj;
						
					default:
						return null;				
				}
				
			case CHAR:
				switch (targetType) {
					case BOOLEAN:
						return ((char) obj) == '1' ? true : false; 
						
					case BYTE:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case CHAR:
						return (char) obj;
						
					case DOUBLE:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
												
					case FLOAT:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
												
					case INT:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case LONG:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case SHORT:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case STRING:
						return obj.toString();
						
					case OBJECT:
						return obj;
						
					default:
						return null;				
				}
					
			case DOUBLE:
				switch (targetType) {
					case BOOLEAN:
						return ((double) obj) == 1.0 ? true : false;
						
					case BYTE:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case CHAR:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case DOUBLE:
						return (double) obj;
						
					case FLOAT:
						return ((Double) obj).floatValue();
						
					case INT:
						return ((Double) obj).intValue();
						
					case LONG:
						return ((Double) obj).longValue();
						
					case SHORT:
						return ((Double) obj).shortValue();
						
					case STRING:
						return obj.toString();
						
					case OBJECT:
						return obj;
						
					default:
						return null;				
				}
				
			case FLOAT:
				switch (targetType) {
					case BOOLEAN:
						return ((float) obj) == 1.0f ? true : false;
						
					case BYTE:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case CHAR:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case DOUBLE:
						return (double) ((float) obj);
						
					case FLOAT:
						return (float) obj;
						
					case INT:
						return ((Float) obj).intValue();
						
					case LONG:
						return ((Float) obj).longValue();
						
					case SHORT:
						return ((Float) obj).shortValue();
						
					case STRING:
						return obj.toString();
						
					case OBJECT:
						return obj;
						
					default:
						return null;				
				}
				
			case INT:
				switch (targetType) {
					case BOOLEAN:
						return ((int) obj) == 1 ? true : false;
						
					case BYTE:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
												
					case CHAR:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case DOUBLE:
						return (double) ((int) obj);
						
					case FLOAT:
						return (float) ((int) obj);
						
					case INT:
						return (int) obj;
						
					case LONG:
						return (long) ((int) obj);
						
					case SHORT:
						return (short) ((int) obj);
						
					case STRING:
						return obj.toString();
						
					case OBJECT:
						return obj;
						
					default:
						return null;				
				}
				
			case LONG:
				switch (targetType) {
					case BOOLEAN:
						return ((long) obj) == 1L ? true : false;
						
					case BYTE:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case CHAR:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case DOUBLE:
						return (double) ((long) obj);
						
					case FLOAT:
						return (float) ((long) obj);
						
					case INT:
						return (int) ((long) obj);
						
					case LONG:
						return (long) obj;
						
					case SHORT:
						return (short) ((long) obj);
						
					case STRING:
						return obj.toString();
						
					case OBJECT:
						return obj;
						
					default:
						return null;				
				}
				
			case SHORT:
				switch (targetType) {
					case BOOLEAN:
						return ((short) obj) == 1 ? true : false;
						
					case BYTE:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case CHAR:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case DOUBLE:
						return (double) ((short) obj);
						
					case FLOAT:
						return (float) ((short) obj);
						
					case INT:
						return (int) ((short) obj);
						
					case LONG:
						return (long) ((short) obj);
						
					case SHORT:
						return (short) obj;
						
					case STRING:
						return obj.toString();
						
					case OBJECT:
						return obj;
						
					default:
						return null;				
				}
				
			case STRING:
				switch (targetType) {
					case BOOLEAN:
						return Boolean.parseBoolean((String) obj);
						
					case BYTE:
						return Byte.parseByte((String) obj);
						
					case CHAR:
						throw new IllegalArgumentException("Conversion from " + sourceType + " to " + targetType + " is not implemented");
						
					case DOUBLE:
						// convert , to .
						s = (String) obj;
						if (s.contains(",")) {
							s = s.replace(".", "");
							s = s.replace(",", ".");
						}
						return Double.parseDouble(s);
						
					case FLOAT:
						// convert , to .
						s = (String) obj;
						s = (String) obj;
						if (s.contains(",")) {
							s = s.replace(".", "");
							s = s.replace(",", ".");
						}
						return Float.parseFloat(s);
						
					case INT:
						// convert , to . TODO this is not safe and does not cover all cases
						s = (String) obj;
						s = s.replace(",", ".");
						// cut string of decimals
						if (s.contains(".")) {
							s = s.substring(0, s.indexOf("."));
							return Integer.parseInt(s);
						} else {
							return Integer.parseInt(s);
						}						
						/* ALTERNATIVE CODE: seems to be slower by a couple ns	
						s = (String) obj;
						Double d = Double.parseDouble(s.replace(",", "."));
						return d.intValue();
						*/						
						
					case LONG:
						// convert , to . TODO this is not safe and does not cover all cases
						s = (String) obj;
						s = s.replace(",", ".");
						// cut string of decimals
						if (s.contains(".")) {
							s = s.substring(0, s.indexOf("."));
							return Long.parseLong(s);
						} else {
							return Long.parseLong(s);
						}
						
					case SHORT:
						// convert , to . TODO this is not safe and does not cover all cases
						s = (String) obj;
						s = s.replace(",", ".");
						// cut string of decimals
						if (s.contains(".")) {
							s = s.substring(0, s.indexOf("."));
							return Short.parseShort(s);
						} else {
							return Short.parseShort(s);
						}
						
					case STRING:
						return obj.toString();
						
					case OBJECT:
						return obj;
						
					default:
						return null;				
				}
			
			case OBJECT:
				switch (targetType) {
					case BOOLEAN:
						return (boolean) obj;
						
					case BYTE:
						return (byte[]) obj;
						
					case CHAR:
						return (char[]) obj;
						
					case DOUBLE:
						return (double) obj;
						
					case FLOAT:
						return (float) obj;
						
					case INT:
						return (int) obj;
						
					case LONG:
						return (long) obj;
						
					case SHORT:
						return (short) obj;
						
					case STRING:
						return obj.toString();
						
					case OBJECT:
						return obj;
						
					default:
						return null;				
				}
				
			default:
				return null;		
		}
	}
	
	/**
	 * converts the objects in {@code objs} from assumed {@code sourceType} to {@code targetType}
	 * TODO this code is not optimized for speed, as convertToDataType is run through for every object in {@code objs}
	 * @param obj
	 * @param sourceType
	 * @param targetType
	 * @return
	 */
	public static Object[] cast(Object[] objs, DataType sourceType, DataType targetType) {
		Object[] newObjs = new Object[objs.length];
		int i = 0;
		for (Object o : objs) {
			newObjs[i] = cast(o, sourceType, targetType);
			++i;
		}
		return newObjs;
	}
	
	/**
	 * converts the objects in {@code objs} from assumed {@code sourceType} to {@code targetType}
	 * TODO this code is not optimized for speed, as convertToDataType is run through for every object in {@code objs}
	 * @param obj
	 * @param sourceType
	 * @param targetType
	 * @return
	 */
	public static Object[] cast(Object[] objs, DataType targetType) {
		Object[] newObjs = new Object[objs.length];
		int i = 0;
		for (Object o : objs) {
			newObjs[i] = cast(o, targetType);
			++i;
		}
		return newObjs;
	}
	
	/**
	 * converts the objects in {@code objs} from assumed {@code sourceType} to {@code targetType}
	 * TODO this code is not optimized for speed, as convertToDataType is run through for every object in {@code objs}
	 * @param objs
	 * @param sourceType
	 * @param targetType
	 * @return
	 */
	public static List<Object> cast(List<Object> objs, DataType sourceType, DataType targetType) {
		List<Object> newObjs = new ArrayList<Object>(objs.size());
		for (Object o : objs) {
			newObjs.add(cast(o, sourceType, targetType));
		}
		return newObjs;
	}
	
	/**
	 * converts the {@code Objects[]} in {@code objs} to double[] assuming {@code DataType} {@code sourceType} 
	 * @param objs
	 * @param sourceType
	 * @return
	 */
	public static double[] castToDouble(Object[] objs, DataType sourceType) {
		double[] ar = new double[objs.length];
		for (int i = 0; i < objs.length; i++) {
			ar[i] = (double) cast(objs[i], sourceType, DataType.DOUBLE);
		}
		return ar;		
	}
	
	/**
	 * converts the {@code Objects[]} in {@code objs} to double[] respect to the underlying {@code objs'} classes
	 * @param objs
	 * @return
	 */
	public static double[] castToDouble(Object[] objs) {
		double[] ar = new double[objs.length];
		for (int i = 0; i < objs.length; i++) {
			ar[i] = (double) cast(objs[i], of(objs[i]), DataType.DOUBLE);
		}
		return ar;
	}
	
	
	/**
	 * returns the (assumed) {@code DataType} of given String {@code value}
	 * @param value
	 * @return
	 */
	public static DataType of(String value) {
		if (value == null) {
			return null;
		} else if (isBoolean(value)) {
			return DataType.BOOLEAN;
		} else if (isInteger(value)) {
			return DataType.INT;
		} else if (isLong(value)) {
			return DataType.LONG;
		} else if (isDouble(value)) {
			return DataType.DOUBLE;
		} else {
			return DataType.OBJECT;
		}
	}
	
	/**
	 * returns the {@code DataType} of this {@code obj}'s {@code Class}
	 * @param obj
	 * @return
	 */
	public static DataType of(Object obj) {
		return of(obj.getClass());
	}

	/**
	 * methods attempts to parse {@code obj} into {@code DataType}, attempting brute force casting
	 * <br>casting is attempted in the following order:
	 * <br>boolean -> integer -> long -> double -> String --> Object
	 * <br>if the {@code obj} cannot be parsed into a primitive, than it is returned as {@code String} or finally {@code Object}
	 * @param str
	 * @return
	 */
	public static Object cast(Object obj) {
		if (obj == null) {
			return null;
		}
		if (obj instanceof String) {
			String str = (String) obj;
			if (isBoolean(str)) {
				return parseBoolean(str);
			}
			if (isInteger(str)) {
				return Integer.parseInt(str);
			}
			if (isLong(str)) {
				return Long.parseLong(str);
			}
			if (isDouble(str)) {
				return Double.parseDouble(str);
			}
			return str;	
		} else {
			return obj;
		}	
	}
	
	/**
	 * return true/false if the specified string can be parsed as {@code boolean}
	 * @param str
	 * @return
	 */
	public static boolean isBoolean(String str) { 
		if (str == null) {
	        return false;
	    }
	    if (str.toLowerCase().contentEquals("yes") || str.toLowerCase().contentEquals("true") || str.toLowerCase().contentEquals("ja") || str.toLowerCase().contentEquals("1") || str.toLowerCase().contentEquals("y") || str.toLowerCase().contentEquals("j") || str.toLowerCase().contentEquals("ok") || str.toLowerCase().contentEquals("i.o.")) {
	    	return true;
	    } else if (str.toLowerCase().contentEquals("no") || str.toLowerCase().contentEquals("false") || str.toLowerCase().contentEquals("nein") || str.toLowerCase().contentEquals("0") || str.toLowerCase().contentEquals("n") || str.toLowerCase().contentEquals("n.i.o.")) {
	    	return true;
	    } else {
	    	return false;
	    }
	}
	
	/**
	 * return true/false if the specified string can be parsed into a double
	 * @param str
	 * @return
	 */
	public static boolean isDouble(String str) { 
		if (str == null) {
	        return false;
	    }
	    try {
	    	if (str.contains(",")) {
	    		str = str.replace(".", "");
	    		str = str.replace(",", ".");
			}
	        double d = Double.parseDouble(str);
	    } catch (NumberFormatException nfe) {
	        return false;
	    }
	    return true;
	}
	
	/**
	 * return true/false if the specified string can be parsed into a float
	 * @param str
	 * @return
	 */
	public static boolean isFloat(String str) { 
		if (str == null) {
	        return false;
	    }
	    try {
	    	if (str.contains(",")) {
	    		str = str.replace(".", "");
	    		str = str.replace(",", ".");
			}
	        double f = Float.parseFloat(str);
	    } catch (NumberFormatException nfe) {
	        return false;
	    }
	    return true;
	}
	
	
	/**
	 * return true/false if the specified string can be parsed into a {@code Integer}
	 * @param str
	 * @return
	 */
	public static boolean isInteger(String str) {
		if (str == null) {
	        return false;
	    }
	    try {
	        int i = Integer.parseInt(str);
	    } catch (NumberFormatException nfe) {
	        return false;
	    }
	    return true;
	}
	
	/**
	 * return true/false if the specified string can be parsed into a {@code Long}
	 * @param str
	 * @return
	 */
	public static boolean isLong(String str) {
		if (str == null) {
	        return false;
	    }
	    try {
	        long lo = Long.parseLong(str);
	    } catch (NumberFormatException nfe) {
	        return false;
	    }
	    return true;
	}
			
	/**
	 * parses the given {@code str} as {@code boolean}
	 * <br>if {@code str} cannot be parsed an {@code IllegalArgumentException} is thrown
	 * @param str
	 * @return
	 */
	private static boolean parseBoolean(String str) {
		if (str.toLowerCase().contentEquals("yes") || str.toLowerCase().contentEquals("true") || str.toLowerCase().contentEquals("ja") || str.toLowerCase().contentEquals("1") || str.toLowerCase().contentEquals("y") || str.toLowerCase().contentEquals("j") || str.toLowerCase().contentEquals("ok") || str.toLowerCase().contentEquals("i.o.")) {
	    	return true;
	    } else if (str.toLowerCase().contentEquals("no") || str.toLowerCase().contentEquals("false") || str.toLowerCase().contentEquals("nein") || str.toLowerCase().contentEquals("0") || str.toLowerCase().contentEquals("n") || str.toLowerCase().contentEquals("n.i.o.")) {
	    	return false;
	    } else {
	    	throw new IllegalArgumentException(str + " cannot be parsed as " + Boolean.class.getSimpleName());
	    }
	}
		
	/**
	 * casts the value type of {@code map} to the given {@code newType}
	 * @param <K, V>
	 * @param map
	 * @param newType
	 * @return
	 */
	public static <K, V> Map castMapValueType(Map<K, V> map, DataType newType) {
		Map newMap;
		switch (newType) {
			case BOOLEAN:
				newMap = new HashMap<String, Boolean>();
				break;
				
			case BYTE:
				newMap = new HashMap<String, Byte>(); 
				break;
				
			case CHAR:
				newMap = new HashMap<String, Character>();
				break;
				
			case DOUBLE:
				newMap = new HashMap<String, Double>();
				break;
				
			case FLOAT:
				newMap = new HashMap<String, Float>();
				break;
				
			case INT:
				newMap = new HashMap<String, Integer>();
				break;
				
			case LONG:
				newMap = new HashMap<String, Long>(); 
				break;
				
			case OBJECT:
				newMap = new HashMap<String, Object>();
				break;
				
			case SHORT:
				newMap = new HashMap<String, Short>();
				break;
				
			case STRING:
				newMap = new HashMap<String, String>();
				break;
				
			default:
				newMap = null;		
		}
		if (newMap != null) {
			for (Entry<K, V> entry : map.entrySet()) {
				newMap.put(entry.getKey(), cast(entry.getValue(), of(entry.getValue()), newType));
			}			
		}
		return newMap;
	}
	
	/**
	 * returns the {@code DataType} of this {@code clazz}
	 * @param clazz
	 * @return
	 */
	public static DataType of(Class<?> clazz) {
		if (clazz.equals(Double.class) || clazz.equals(double.class) || clazz.equals(Double[].class) || clazz.equals(double[].class)) {
			return DOUBLE;
		} else if (clazz.equals(Integer.class) || clazz.equals(int.class) || clazz.equals(Integer[].class) || clazz.equals(int[].class)) {
			return INT;
		} else if (clazz.equals(Short.class) || clazz.equals(short.class) || clazz.equals(Short[].class) || clazz.equals(short[].class)) {
			return SHORT;
		} else if (clazz.equals(Long.class) || clazz.equals(long.class) || clazz.equals(Long[].class) || clazz.equals(long[].class)) {
			return LONG;
		} else if (clazz.equals(Character.class) || clazz.equals(char.class) || clazz.equals(Character[].class) || clazz.equals(char[].class)) {
			return CHAR;
		} else if (clazz.equals(String.class) || clazz.equals(String[].class)) {
			return STRING;
		} else if (clazz.equals(Float.class) || clazz.equals(float.class) || clazz.equals(Float[].class) || clazz.equals(float[].class)) {
			return FLOAT;
		} else if (clazz.equals(Boolean.class) || clazz.equals(boolean.class) || clazz.equals(Boolean[].class) || clazz.equals(boolean[].class)) {
			return BOOLEAN;
		} else {
			return OBJECT;
		}
	}
	
	public static DataType getIntTypeByBytes(int bytes) {
		switch (bytes) {
			case 1:
				return DataType.BYTE;
			
			case 2:
				return DataType.SHORT;
				
			case 4:
				return DataType.INT;
				
			case 8:
				return DataType.LONG;
						
			default:
				return null;
		}		
	}
			
	public void setBytes(int bytes) {
		this.bytes = bytes;
	}
		
	@Override
	public String toString() {
		return this.toString(true);
	}
	
	public String toString(boolean showBytes) {
		if (showBytes) {
			return super.toString() + "(" + this.bytes + ")";
		} else {
			return super.toString();
		}
	}
	
	public int nBytes() {
		return this.bytes;
	}

	public static DataType CHAR(int i) {
		DataType dt = DataType.CHAR;
		dt.setBytes(i);
		return dt;
	}
	
	public static DataType BYTE(int i) {
		DataType dt = DataType.BYTE;
		dt.setBytes(i);
		return dt;
	}
	
	public static DataType STRING(int i) {
		DataType dt = DataType.STRING;
		dt.setBytes(i);
		return dt;
	}
	
	public static Map<String, DataType> toMap(){
		Map<String, DataType> map = new LinkedHashMap<String, DataType>();
		for (DataType dt : DataType.values()) {
			map.put(dt.toString(), dt);
		}
		return map;
	}
}
