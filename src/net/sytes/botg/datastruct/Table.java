package net.sytes.botg.datastruct;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;

import net.sytes.botg.array.Ar;
import net.sytes.botg.datatypes.DataType;

/**
 * Table implementation in Java based on a {@code LinkedHashMap} for addressing columns by name, the data within a column is represented by {@code ArrayList<Object>} for easy adding and removing of rows
 * <br>indices of row and column start with 0 
 * <br>the implementation provides functionality like xlookup, removing duplicates or filter functions known from excel tables
 */
public class Table implements Cloneable {
	
	/**
	 * data container for the columnwise data
	 */
	private Map<String, List<Object>> data = new LinkedHashMap<String, List<Object>>();	
	
	/**
	 * defines cloningBehavior of data coming and out of this {@code Table} 
	 */
	private CloningBehavior cloningBehavior;
	
	/**
	 * default cloning behavior config
	 */
	private static final CloningBehavior DEFAULT_CLONING_BEHAVIOR = CloningBehavior.CLONE_ON_ENTRY; 
	
	public enum CloningBehavior {
		CLONE_ON_ENTRY, ALWAYS_CLONE, NONE
	}
	
	public enum Comparator {
		
		EQUAL("="), SMALLER_THAN("<"), GREATER_THAN(">"), EQUAL_OR_SMALLER("<="), EQUAL_OR_GREATER(">="), LIKE("~"), NOT_NULL("!=NULL"), NOT_EQUAL("!=");
		
		private String comperator;
		
		private Comparator(String comperator) {
			this.comperator = comperator;
		}
		
		public String toString() {
			return this.comperator;
		}
		
		public static Comparator containsComperator(String s) {
		    if (s.contains(NOT_NULL.toString())) {
		    	return NOT_NULL;
			} else if (s.contains(NOT_EQUAL.toString())) {
				return NOT_EQUAL;
			} else if (s.contains(EQUAL_OR_SMALLER.toString())) {
				return EQUAL_OR_SMALLER;
			} else if (s.contains(EQUAL_OR_GREATER.toString())) {
				return EQUAL_OR_GREATER;
			} else if (s.contains(EQUAL.toString())) {
				return EQUAL;
			} else if (s.contains(SMALLER_THAN.toString())) {
				return SMALLER_THAN;
			} else if(s.contains(GREATER_THAN.toString())) {
				return GREATER_THAN;
			} else if (s.contains(LIKE.toString())) {
				return LIKE;
			} else {
				return null;
			}
		}
		
		public static Comparator toEnum(String comparator) {
			switch (comparator) {
				case "=":
					return EQUAL;
				
				case "!=":
					return NOT_EQUAL;
					
				case "<":
					return SMALLER_THAN;
					
				case "<=":
					return EQUAL_OR_SMALLER;
					
				case ">":
					return GREATER_THAN;
					
				case ">=":
					return EQUAL_OR_GREATER;
				
				case "~":
					return LIKE;
					
				case "!=NULL":
					return NOT_NULL;
				
				default:
					return null;
			}
		}	
		
	}
	
	public enum MatchMode {
		EXACT_MATCH, EXACT_OR_NEXT_SMALLER, EXACT_OR_NEXT_LARGER, WILDCARD_MATCH
	}
	
	public enum SearchMode {
		START_WITH_FIRST, START_WITH_LAST, ASCENDING, DESCENDING
	}
	
	public enum SortMode {
		ASCENDING, DESCENDING, A_TO_Z, Z_TO_A
	}
	
	public Table() {
		this(null, null, DEFAULT_CLONING_BEHAVIOR);
	}
	
	public Table(CloningBehavior cloningBehavior) {
		this(null, null, cloningBehavior);
	}
	
	public Table(String columnName, List<Object> objs, CloningBehavior cloningBehavior) {
		this.cloningBehavior = cloningBehavior;
		this.addColumn(columnName, objs);
	}
	
	/**
	 * the object array is added as a column with default name col1
	 * @param ar
	 */
	public Table(Object[] ar) {
		this("COL0", Arrays.asList(ar), DEFAULT_CLONING_BEHAVIOR);
	}
	
	/**
	 * the object array is added as a column with columnName
	 * @param columnName
	 * @param ar
	 */
	public Table(String columnName, Object[] ar) {
		this(columnName, Arrays.asList(ar), DEFAULT_CLONING_BEHAVIOR);
	}
		
	/**
	 * the list is stored as columnName in data
	 * @param columnName
	 * @param data
	 */
	public Table(String columnName, List<Object> data) {
		this(columnName, data, DEFAULT_CLONING_BEHAVIOR);
	}
		
	/**
	 * the map is interpreted as Columnwise data, where the keys specify the name of the column and the values are Objects or List<Object>'s containing the row data
	 * @param mapData
	 */	
	public Table(Map<String, Object> mapData) {
		this(null, null, DEFAULT_CLONING_BEHAVIOR);
		for (Entry<String, Object> entry : mapData.entrySet()) {
			if(Ar.isArray(entry.getValue())) {
				this.addColumn(entry.getKey(), (List<Object>) entry.getValue());
			} else if (entry.getValue() instanceof List<?>) {
				this.addColumn(entry.getKey(), (List<Object>) entry.getValue());
			} else {
				this.add(entry.getKey(), entry.getValue());
			}
		}
	}
	
	/**
	 * the array is interpreted as table in the following way for number of columns and rows: 
	 * <br>ar[numOfColumns][numOfRows], length of header String and number of columns must match
	 * <br>the columns are named Col0, Col1, Col2, ...
	 * @param ar
	 */
	public Table(Object[][] ar) {		
		this(DEFAULT_CLONING_BEHAVIOR);
		for (int i = 0; i < ar.length; i++) {
			this.addColumn("COL" + (i), Arrays.asList(ar[i]));
		}
	}
	
	/**
	 * the array is interpreted as table in the following way for number of columns and rows: 
	 * <br>ar[numOfColumns][numOfRows], length of header String and number of columns must match
	 * <br>strings in the header are used for columnNames
	 * 
	 * @param headers
	 * @param ar
	 */
	public Table(String[] headers, Object[][] ar) {
		this(DEFAULT_CLONING_BEHAVIOR);
		if (headers.length != ar.length) {
			throw new IllegalArgumentException("length of headers array (" + headers.length + ") and number of columns in ar (" + ar.length + ") must match");
		}
		int i = 0;
		for (String h : headers) {
			this.addColumn(h, Arrays.asList(ar[i]));
			++i;
		}
	}

	/**
	 * the array is interpreted as table in the following way for number of columns and rows: 
	 * <br>ar[numOfColumns][numOfRows], length of header String and number of columns must match 
	 * @param headers List of Strings
	 * @param ar 2D array containing columns with data
	 */
	public Table(List<String> headers, Object[][] ar) {
		this(DEFAULT_CLONING_BEHAVIOR);
		if (headers.size() != ar.length) {
			throw new IllegalArgumentException("length of headers list (" + headers.size() + ") and number of columns in ar (" + ar.length + ") must match");
		}
		int i = 0;
		for (String h : headers) {
			this.data.put(h, Arrays.asList(ar[i]));
			++i;
		}
	}
	
	/**
	 * initializes the Table with empty Columns with {@code headers} specified
	 * @param headers
	 */
	public Table(String[] headers){
		this(DEFAULT_CLONING_BEHAVIOR);
		for (String h : headers) {
			this.data.put(h, new ArrayList<Object>());
		}
	}
	
	/**
	 * initializes the Table with empty Columns with the {@code headers} specified
	 * @param headers
	 */
	public Table(List<String> headers){
		this(DEFAULT_CLONING_BEHAVIOR);
		for (String h : headers) {
			this.data.put(h, new ArrayList<Object>());
		}
	}
	
	/**
	 * returns the value at {@code column} at row {@code r}
	 * @param column
	 * @param r
	 * @return
	 */
	public Object get(String column, int r) {
		this.checkColumn(column);
		this.checkRange(r);
		return this.getColumn(column).get(r);
	}
	
	/**
	 * returns the value at row {@code r} in column {@code c}
	 * @param r
	 * @param c
	 * @return
	 */
	public Object get(int r, int c) {
		this.checkRange(r, c);
		return this.getColumn(c).get(r);
	}
	

	/**
	 * returns a new sub {@code Table} containing all the columns specified by {@code columnNames}
	 * @param columnNames
	 * @return
	 */
	public Table get(String[] columnNames) {
		Table t = new Table();	
		for (String columnName : columnNames) {
			// check if columnName is present in this table
			this.checkRange(0, 0, columnName);
			t.addColumn(columnName, this.getColumn(columnName));
		}	
		return t;
	}
	
	/**
	 * return row with index {@code i} as new {@code Table}
	 * @param i
	 * @return
	 */
	public Table getRowAsTable(int i) {
		Table t = new Table();
		Map<String, Object> map = this.getRowMap(i);
		t.addRow(map);
		return t;
	}
	
	/**
	 * returns row values specified by {@code columnNames} and {@code row} index
	 * @param columnNames
	 * @param row
	 * @return
	 */
	public List<Object> get(String[] columnNames, int row) {
		return this.get(Arrays.asList(columnNames), row);	
	}
	
	/**
	 * returns row values specified by {@code columnNames} and {@code row} index
	 * @param columnNames
	 * @param row
	 * @return
	 */
	public List<Object> get(List<String> columnNames, int row) {
		this.checkColumn(columnNames);
		this.checkRange(row);
		List<Object> rowData = new ArrayList<Object>();
		for (String columnName : columnNames) {
			rowData.add(this.data.get(columnName).get(row));
		}
		return rowData;		
	}

	/**
	 * returns a new sub {@code Table} containing all the columns specified by {@code columnIndices}
	 * @param columns
	 * @return
	 */
	public Table get(int[] columnIndices) {
		// check if all column indices are in range
		Table t = new Table();	
		for (int c : columnIndices) {
			this.checkRange(0, c);
			String columnName = this.getColumnName(c);
			t.addColumn(columnName, this.getColumn(columnName));
		}		
		return t;
	}
	
	/**
	 * get the row data as list of objects specified by r
	 * @param r
	 * @return
	 */
	public List<Object> getRow(int r){
		this.checkRange(r);
		List<Object> rowData = new ArrayList<Object>();
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			rowData.add(entry.getValue().get(r));
		}
		return rowData;
	}
	
	/**
	 * get the specified row {@code r} data as a {@code Map<String, Object>}, where the key holds the column names
	 * @param r
	 * @return
	 */
	public Map<String, Object> getRowMap(int r){
		this.checkRange(r);
		Map<String, Object> row = new LinkedHashMap<String, Object>();
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			row.put(entry.getKey(), entry.getValue().get(r));
		}
		return row;
	}
	
	/**
	 * get the column data as list of objects specified by c
	 * <br>this creates a deep clone of the underlying objects inside the list
	 * @param c
	 * @return
	 */
	public List<Object> getColumn(int c){
		int i = 0;
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			if (c == i) {
				return this.getColumn(entry.getKey());
			}
			++i;
		}
		return null;
	}
	
	/**
	 * get the column data as list of objects specified by {@code columnName}
	 * <br>this creates a clone of the {@code ArrayList} containing this column's elements
	 * @param columnName
	 * @return
	 */
	public List<Object> getColumn(String columnName){
		List<Object> newList = new ArrayList<Object>(this.data.get(columnName).size());
		//Collections.copy(newList, this.data.get(columnName));
		newList.addAll(this.data.get(columnName));
		return newList;
	}
	
	/**
	 * return all columNames
	 * @return
	 */
	public String[] getColumnNames() {
		String[] columnNames = this.data.keySet().toArray(new String[] {});
		return columnNames;
	}
	
	/**
	 * return columnName at column index {@code c}
	 * @param c
	 * @return
	 */
	public String getColumnName(int c) {
		this.checkRange(0, c);
		String columnName = this.data.keySet().toArray(new String[] {})[c];
		return columnName;
	}
	
	/**
	 * add empty column with columnName
	 * @param columnName
	 */
	public void addColumn(String columnName) {
		if (this.getNumberOfElementsInRows() == 0) {
			this.data.put(columnName, new ArrayList<Object>());
		} else {
			this.addColumn(columnName, new ArrayList<Object>());
		}
	}

	/**
	 * add column with columnData
	 * @param columnData
	 */
	public void addColumn(Object[] columnData) {
		this.addColumn("COL" + (this.columns()), columnData);
	}
	
	/**
	 * add a new column with columnName with columnData
	 * @param columnName
	 * @param columnData
	 */
	public void addColumn(String columnName, Object[] columnData) {
		List<Object> newList = new ArrayList<Object>();
		newList.addAll(Arrays.asList(columnData));
		this.addColumn(columnName, newList);
	}
	
	public void addColumn(List<Object> columnData) {
		String columnName = "COL" + this.columns();
		this.addColumn(columnName, columnData);
	}
	
	public void addColumn(String columnName, List<Object> columnData) {
		if (columnName == null) {
			return;
		}
		if (columnData == null) {
			return;
		}
		if (this.hasColumn(columnName)) {
			throw new IllegalArgumentException("a column with name '" + columnName + "' is already present. Try replace or rename the column you want to add.");
		}
		List<Object> clonedData;
		switch (this.cloningBehavior) {
			case ALWAYS_CLONE:
				clonedData = new ArrayList<Object>(columnData.size());
				//Collections.copy(clonedData, columnData);
				clonedData.addAll(columnData);
				this.data.put(columnName, clonedData);
				break;
				
			case CLONE_ON_ENTRY:
				clonedData = new ArrayList<Object>(columnData.size());
				//Collections.copy(clonedData, columnData);
				clonedData.addAll(columnData);
				this.data.put(columnName, clonedData);
				break;
				
			case NONE:
				this.data.put(columnName, columnData);
				break;
				
			default:
				throw new IllegalStateException("");		
		}
		this.fillUpColumns();
	}
	
	/**
	 * adds a new row based on the array rowData, if rowData does not equal available columns an Exception is thrown 
	 * @param rowData
	 */
	public void addRow(Object[] rowData) {
		List<Object> newList = new ArrayList<Object>();
		newList.addAll(Arrays.asList(rowData));
		this.addRow(newList);
	}
	
	/**
	 * adds a new row based on the data in rowData, if rowData does not equal available columns an Exception is thrown
	 * @param rowData
	 */
	public void addRow(List<Object> rowData) {
		if (rowData.size() != this.columns()) {
			throw new IllegalArgumentException("elements in rowData (" + rowData.size() + ") does not equal the number of columns (" + this.columns() + ") in this Table");
		}
		/*
		if (this.cloneOnEntry) {
			List<Object> clonedData = new ArrayList<Object>(rowData.size());
			//Collections.copy(clonedData, rowData);
			clonedData.addAll(rowData);
			int c = 0;
			for (Entry<String, List<Object>> entry : this.data.entrySet()) {
				entry.getValue().add(clonedData.get(c));
				++c;
			}
		} else {
			int c = 0;
			for (Entry<String, List<Object>> entry : this.data.entrySet()) {
				entry.getValue().add(rowData.get(c));
				++c;
			}
		}
		*/
		
		// TODO here should be some switch case with different cloning behavior
		int c = 0;
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			entry.getValue().add(rowData.get(c));
			++c;
		}
	}
	
	/**
	 * adds elements within a new row based on data in Map
	 * @param rowData
	 */
	public void addRow(Map<String, Object> rowData) {
		for (Entry<String, Object> entry : rowData.entrySet()) {			
			this.data.get(entry.getKey()).add(entry.getValue());
		}
		this.fillUpColumns();
	}
	
	/**
	 * adds a new element in columnName
	 * @param columnName
	 * @param value
	 */
	public void add(String columnName, Object value) {
		if (!data.containsKey(columnName)) {
			this.data.put(columnName, new ArrayList<Object>());
		}
		this.data.get(columnName).add(value);
		this.fillUpColumns();
	}
	
	/**
	 * adds a new {@code value} to column specified by index {@code col}
	 * @param col
	 * @param value
	 */
	public void add(int col, Object value) {
		int c = 0;
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			if (c == col) {
				entry.getValue().add(value);
			}
			++c;
		}
		this.fillUpColumns();
	}
	
	/**
	 * adds all columns from {@code table} to this {@code Table} if overwrite is set to true,
	 * <br>then the column with the same name in {@code table} is kept
	 * @param table
	 */
	public void addColumns(Table table, boolean overwrite) {
		for (Entry<String, List<Object>> entry : table.data.entrySet()) {
			this.data.put(entry.getKey(), entry.getValue());
		}
		this.fillUpColumns();
	}
	
	/**
	 * adds all rows of table to the bottom of this object, if the number of columns is equal
	 * @param table
	 */
	public void addRows(Table table) {
		this.addRows(table, false);
	}
		
	/**
	 * Adds all rows of table to the bottom of this object, if the number of columns is equal.
	 * <br>If ignoreColumnNames is false, then only rows of columns from table present in this object are added,
	 * <br>else the rows are added to the same index of columns they came from in table
	 * @param table
	 * @param ignoreColumnNames
	 */
	public void addRows(Table table, boolean ignoreColumnNames) {
		if (ignoreColumnNames) {
			int i = 0;
			for (Entry<String, List<Object>> entry : this.data.entrySet()) {
				int j = 0;
				for (Entry<String, List<Object>> entry2 : table.data.entrySet()) {
					if (i == j) {
						entry.getValue().addAll(entry2.getValue());
					}
					++j;
				}
				++i;
			}
		} else {
			for (Entry<String, List<Object>> entry : table.data.entrySet()) {
				this.data.get(entry.getKey()).addAll(entry.getValue());
			}
		}
		this.fillUpColumns();
	}
	
	/**
	 * replace an element at column specified by {@code columnName} and {@code row} index with {@code value}
	 * @param columnName
	 * @param index
	 * @param value
	 */
	public void set(String columnName, int row, Object value) {
		this.checkRange(row, 0, columnName);
		this.data.get(columnName).set(row, value);
	}
	
	/**
	 * replace the whole row with new data
	 * @param row
	 * @param data
	 */
	public void set(int row, Object[] data) {
		if (data.length != this.columns()) {
			throw new IllegalArgumentException("elements in rowData (" + data.length + ") does not equal the number of columns (" + this.columns() + ") in this Table");
		}
		int c = 0;
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			entry.getValue().set(row, data[c]);
			++c;
		}
	}

	/**
	 * replace an element specified by (row, column) with value
	 * @param row
	 * @param column
	 * @param value
	 */
	public void set(int row, int column, Object value) {
		String columnName = this.getColumnName(column);
		this.set(columnName, row, value);
	}
	
	/**
	 * removes the value specified with {@code columnName} and {@code row} index
	 * @param columnName
	 * @param row
	 */
	public void remove(String columnName, int row) {
		// check if indices are within limits
		this.checkRange(row, 0, columnName);
		this.nullifyElement(columnName, row);
		if (this.isRowRemovable(row)) {
			this.removeRow(row);
		}
	}
	
	/**
	 * removes the whole row at {@code row} index
	 * @param row
	 */
	public void remove(int row) {
		// check if indices are within limits
		this.checkRange(row);
		this.removeRow(row);
	}
	
	/**
	 * removes the value specified at (row, column) in Table
	 * @param row
	 * @param column
	 */
	public void remove(int row, int column) {
		// check if indices are within limits
		this.checkRange(row, column);
		this.nullifyElement(row, column);
		if (this.isRowRemovable(row)) {
			this.removeRow(row);
		}
	}
	
	/**
	 * removes a column specified by {@code columnName}
	 * @param columnName
	 */
	public void remove(String columnName) {
		this.data.remove(columnName);
	}
	
	/**
	 * removes the column at {@code index}
	 * @param column
	 */
	public void removeColumn(int index) {
		this.checkRange(0, index);
		int i = 0;
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			if (i == index) {
				this.data.remove(entry.getKey());
				return;
			}
			++i;
		}
	}
	
	/**
	 * removes all duplicate rows in this {@code Table} based on all columns in other rows
	 * <br>duplicates are removed from top to bottom, that means entries with lower index are kept
	 */
	public void removeDuplicates() {
		int n = this.getNumberOfElementsInRows();
		int r = 0;
		while (r < n) {
			List<Object> row = this.getRow(r);
			for (int rr = n - 1; rr > r; rr--) {
				List<Object> row2 = this.getRow(rr);
				if (row.equals(row2)) {
					this.remove(rr);
				}
			}
			++r;
			n = this.getNumberOfElementsInRows();
		} 
	}
	
	/**
	 * removes duplicate rows in this {@code Table} based on values in column {@code columnName}
	 * <br>duplicates are removed from top to bottom, that means entries with lower index are kept
	 * @param columnName
	 */
	public void removeDuplicates(String columnName) {
		int n = this.getNumberOfElementsInRows();
		int r = 0;
		List<Object> columnData = this.getColumn(columnName);
		while (r < n) {
			Object value = columnData.get(r);
			for (int rr = n - 1; rr > r; rr--) {
				Object value2 = columnData.get(rr);
				if (value.equals(value2)) {
					this.remove(rr);
				}
			}
			++r;
			n = this.getNumberOfElementsInRows();
		} 
	}
	
	/**
	 * removes duplicate rows in this {@code Table} based on values in columns specified by {@code columnNames}
	 * <br>duplicates are removed from top to bottom, that means entries with lower index are kept
	 * @param columnNames
	 */
	public void removeDuplicates(String[] columnNames) {
		this.removeDuplicates(Arrays.asList(columnNames));
	}
	
	/**
	 * removes duplicate rows in this {@code Table} based on values in columns specified by {@code columnNames}
	 * <br>duplicates are removed from top to bottom, that means entries with lower index are kept
	 * @param columnNames
	 */
	public void removeDuplicates(List<String> columnNames) {
		int n = this.getNumberOfElementsInRows();
		int r = 0;
		while (r < n) {
			List<Object> row = this.get(columnNames, r);
			for (int rr = n - 1; rr > r; rr--) {
				List<Object> row2 = this.get(columnNames,rr);
				if (row.equals(row2)) {
					this.remove(rr);
				}
			}
			++r;
			n = this.getNumberOfElementsInRows();
		} 
	}
			
	public void sort(String sortColumn, SortMode sortMode) {
		// TODO
	}
	
	// TODO
	public Object xLookup(Object searchValue, String searchColumn, String lookupColumn, MatchMode matchMode, SearchMode searchMode) {
		Objects.requireNonNull(searchValue);
		this.checkColumn(searchColumn);
		this.checkColumn(lookupColumn);
		
		List<Object> searchData = new ArrayList<Object>(this.data.get(searchColumn).size());
		//Collections.copy(searchData, this.data.get(searchColumn));
		searchData.addAll(searchData);
		switch (searchMode) {
			case ASCENDING:
				// WARNING: this is not clean --> TODO Fix for unknown object type
				searchData.sort(null);
				break;
				
			case DESCENDING:
				// TODO
				throw new UnsupportedOperationException("SearchMode DESCENDING is not yet implemented");
				
			case START_WITH_FIRST:
				// nothing to do
				break;
				
			case START_WITH_LAST:
				Collections.reverse(searchData);
				break;
				
			default:
				return null;		
		}
		
		int r = 0;
		switch (matchMode) {
			case EXACT_MATCH:				
				for (Object o : searchData) {
					if (o.equals(searchValue)) {
						return this.getColumn(lookupColumn).get(r);
					}
					++r;
				}
				break;
				
			case EXACT_OR_NEXT_LARGER:
				break;
				
			case EXACT_OR_NEXT_SMALLER:
				break;
				
			case WILDCARD_MATCH:
				break;
				
			default:
				return null;		
		}
		return null;
	}
	
	/**
	 * filters this {@code Table}'s column {@code filterColumn} for values that match the filter defined by {@code filterValue} and a specified {@code Comparator comparator}
	 * @param filterValue
	 * @param filterColumn
	 * @param comparator
	 * @return
	 */
	public Table filter(Object filterValue, String filterColumn, Comparator comparator) {
		Table filteredTable = this.clone();
		List<Object> filterData = this.getColumn(filterColumn);
		int n = this.getNumberOfElementsInRows();
		
		String sFilter;
		String sObj;
		
		double fVal;
		double objVal;		
		
		for (int r = n - 1; r >= 0; r--) {
			
			Object obj = filterData.get(r);
			
			switch (comparator) {
				case EQUAL:
					if (!obj.equals(filterValue)) {
						filteredTable.removeRow(r);
					}
					break;
					
				case EQUAL_OR_GREATER:
					fVal = (double) Objects.requireNonNull(DataType.cast(filterValue, DataType.of(filterValue.getClass()), DataType.DOUBLE));				
					objVal = (double) Objects.requireNonNull(DataType.cast(obj, DataType.of(obj.getClass()), DataType.DOUBLE));
					if (objVal < fVal) {
						filteredTable.removeRow(r);
					}
					break;
					
				case EQUAL_OR_SMALLER:
					fVal = (double) Objects.requireNonNull(DataType.cast(filterValue, DataType.of(filterValue.getClass()), DataType.DOUBLE));				
					objVal = (double) Objects.requireNonNull(DataType.cast(obj, DataType.of(obj.getClass()), DataType.DOUBLE));
					if (objVal > fVal) {
						filteredTable.removeRow(r);
					}
					break;
					
				case GREATER_THAN:
					fVal = (double) Objects.requireNonNull(DataType.cast(filterValue, DataType.of(filterValue.getClass()), DataType.DOUBLE));				
					objVal = (double) Objects.requireNonNull(DataType.cast(obj, DataType.of(obj.getClass()), DataType.DOUBLE));
					if (objVal <= fVal) {
						filteredTable.removeRow(r);
					}
					break;
					
				case SMALLER_THAN:
					fVal = (double) Objects.requireNonNull(DataType.cast(filterValue, DataType.of(filterValue.getClass()), DataType.DOUBLE));				
					objVal = (double) Objects.requireNonNull(DataType.cast(obj, DataType.of(obj.getClass()), DataType.DOUBLE));
					if (objVal >= fVal) {
						filteredTable.removeRow(r);
					}
					break;
										
				case LIKE:
					sFilter = (String) Objects.requireNonNull(DataType.cast(filterValue, DataType.of(filterValue.getClass()), DataType.STRING));						
					sObj = (String) Objects.requireNonNull(DataType.cast(obj, DataType.of(obj.getClass()), DataType.STRING));
					if (!sObj.contentEquals(sFilter)) {
						filteredTable.removeRow(r);
					}
					break;
					
				case NOT_EQUAL:
					if (obj.equals(filterValue)) {
						filteredTable.removeRow(r);
					}
					break;
					
				case NOT_NULL:
					if (obj == null) {
						filteredTable.removeRow(r);
					}
					break;
					
				default:
					break;			
			}
		}
		return filteredTable;
	}
	
	/**
	 * nullifies element at specified location given by columnName and row index
	 * @param columnName
	 * @param row
	 */
	private void nullifyElement(String columnName, int row) {
		this.data.get(columnName).set(row, null);
	}
	
	/**
	 * nullifies element at specified location given by (row, column)
	 * @param row
	 * @param column
	 */
	private void nullifyElement(int row, int column) {
		int c = 0;
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			if (c == column) {
				entry.getValue().set(row, null);
				return;
			}
			++c;
		}
	}
	
	/**
	 * checks whether all elements in row are NULL
	 * @param row
	 * @return
	 */
	private boolean isRowRemovable(int row) {
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			if (entry.getValue().get(row) != null) {
				return false;
			}
		}
		return true;
	}
	
	/**
	 * fills up all rows with NULL, if one columns has more rows than the others
	 */
	private void fillUpColumns() { 
		int maxRows = this.getNumberOfElementsInRows();
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			while (entry.getValue().size() < maxRows) {
				entry.getValue().add(null);
			}
		}
	}
	
	/**
	 * returns the number of rows (maximum) in the {@code Table}
	 * @return
	 */
	public int rows() {
		return this.getNumberOfElementsInRows();
	}
	
	/**
	 * returns the largest row number in Table
	 * @return
	 */
	private int getNumberOfElementsInRows() {
		int max = 0;
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			if (entry.getValue().size() > max) {
				max = entry.getValue().size();
			}
		}
		return max;
	}
	
	/**
	 * return the number of columns in Table
	 * @return
	 */
	public int columns() {
		if (this.data.isEmpty()) {
			return 0;
		} else {
			return this.data.values().size();
		}
	}
	
	/**
	 * returns the size of the Table as int array with [rows, cols]
	 * @return
	 */
	public int[] size() {
		int[] s = new int[2];
		s[0] = this.getNumberOfElementsInRows();
		s[1] = this.columns();
		return s;
	}
	
	/**
	 * returns the total number of elements in the table including NULL entries
	 * numOfElems = rows * cols
	 * @return
	 */
	public int getNumberOfElements() {
		return this.getNumberOfElementsInRows() * this.columns();
	}
	
	/**
	 * remove a whole row in Table
	 * @param row
	 */
	private void removeRow(int row) {
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			entry.getValue().remove(row);
		}
	}
		
	/**
	 * returns true/false whether columnName is present in Table
	 * @param columnName
	 * @return
	 */
	public boolean hasColumn(String columnName) {
		if (this.data.containsKey(columnName)) {
			return true;
		} else {
			return false;
		}
	}
		
	/**
	 * checks whether row is within range in Table
	 * @param row
	 */
	public boolean isRowIndexWithinLimit(int row) {
		if (this.getNumberOfElementsInRows() >= row) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * checks whether column is within range in Table
	 * @param column
	 */
	public boolean isColumnIndexWithinLimit(int column) {
		if (this.columns() > column) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * checks whether row index is within Table range
	 * @param row
	 */
	private void checkRange(int row) {
		this.checkRange(row, 0, null);
	}
	
	/**
	 * checks whether row and column index are qithin Table range
	 * @param row
	 * @param col
	 */
	private void checkRange(int row, int col) {
		this.checkRange(row, col, null);
	}
	
	/**
	 * checks whether row, column and columnName are available in Table
	 * @param row
	 * @param col
	 * @param columnName
	 */
	private void checkRange(int row, int col, String columnName) {
		if (!this.isColumnIndexWithinLimit(col)) {
			throw new IndexOutOfBoundsException("column index (" + col + ") is out of range (max index: " + (this.columns() - 1) + ")");
		}
		if (!this.isRowIndexWithinLimit(row)) {
			throw new IndexOutOfBoundsException("row index (" + row + ") is out of range (max index: " + (this.getNumberOfElementsInRows() - 1) + ")");
		}
		if (columnName != null) {
			if (!this.hasColumn(columnName)) {
				throw new IllegalArgumentException("Column (" + columnName + ") is not present in Table!");
			}
		}
	}
	
	/**
	 * checks whether columnName is available in this {@code Table}
	 * <br>if not an {@code IllegalArgumentException} is thrown 
	 * @param columnName
	 */
	private void checkColumn(String columnName) {
		if (!this.data.containsKey(columnName)) {
			throw new IllegalArgumentException(this.getClass().getSimpleName() + " does not contain column '" + columnName + "'.");
		}
	}
	
	/**
	 * checks whether all columnNames are available in this {@code Table}
	 * <br>if not an {@code IllegalArgumentException} is thrown
	 * @param columnNames
	 */
	private void checkColumn(List<String> columnNames) {
		for (String columnName : columnNames) {
			this.checkColumn(columnName);
		}
	}
	
	/**
	 * checks whether all columnNames are available in this {@code Table}
	 * <br>if not an {@code IllegalArgumentException} is thrown
	 * @param columnNames
	 */
	private void checkColumn(String[] columnNames) {
		for (int c = 0; c < columnNames.length; c++) {
			this.checkColumn(columnNames[c]);
		}
	}

	/**
	 * returns the content of this {@code Table} as JSON string
	 * @return
	 */
	public String toJson() {
		StringBuilder sb = new StringBuilder();
		sb.append("{\n");
		int c = 0;
		int cmax = this.data.size();
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {			
			sb.append("\t\"").append(entry.getKey()).append("\": ").append(entry.getValue().toString());
			++c;
			if (c == cmax) {
				sb.append("\n");
			} else {
				sb.append(",\n");
			} 
		}
		sb.append("}");
		return sb.toString();
	}

	/**
	 * returns the table as map
	 * @return
	 */
	public Map<String, List<Object>> toMap() {
		return this.data;
	}
		
	/**
	 * string representation of this Table
	 * @return
	 */
	@Override
	public String toString() {
		String s = Table.class.getSimpleName() + "(" + this.getNumberOfElementsInRows() + ", " + this.columns() + "):\n";
		String temp = this.data.toString();
		temp = temp.replace("], ", "]\n\t");
		temp = temp.replace("{", "{\n\t");
		temp = temp.replace("}", "\n}");
		return s + temp;
	}
	
	/**
	 * on cloning this {@code Table} all columns are added to a new {@code Table}
	 */
	@Override
	public Table clone() {
		Table t = new Table();
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			List<Object> list = new ArrayList<Object>();
			list.addAll(entry.getValue());
			t.addColumn(entry.getKey(), list);
		}
		return t;
	}
	
}

