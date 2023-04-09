package net.sytes.botg.datatypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;

/**
 * Table implementation in Java based on a {@code LinkedHashMap} for adressing columns by name, the data within a column is represented by {@code ArrayList<Object>} for easy adding and removing of rows
 * indices of row and column start with 0 
 * <br>the implementation provides functionality like xlookup, removing duplicates or filter functions known from excel tables
 */
public class Table implements Cloneable {

	private Map<String, List<Object>> data = new LinkedHashMap<String, List<Object>>();	
	
	public enum Comparator {
		EQUALS, SMALLER_THAN, GREATER_THAN, EQUAL_OR_SMALLER, EQUAL_OR_GREATER, LIKE
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
	}
	
	/**
	 * the object array is added as a column with default name col1
	 * @param ar
	 */
	public Table(Object[] ar) {
		this("COL0", ar);
	}
	
	/**
	 * the object array is added as a column with columnName
	 * @param columnName
	 * @param ar
	 */
	public Table(String columnName, Object[] ar) {
		this.data.put(columnName, Arrays.asList(ar));
	}
	
	/**
	 * the object list is stored in data with default columnName
	 * @param data
	 */
	public Table(List<Object> data) {
		this("COL0", data);
	}
	
	/**
	 * the list is stored as columnName in data
	 * @param columnName
	 * @param data
	 */
	public Table(String columnName, List<Object> data) {
		this.data.put(columnName, data);
	}
		
	/**
	 * the map is interpreted as Columnwise data, where the keys specify the name of the column and the values are Object[] arrays containing the row data
	 * @param mapData
	 */
	public Table(Map<String, Object[]> mapData) {
		for (Entry<String, Object[]> entry : mapData.entrySet()) {
			this.data.put(entry.getKey(), Arrays.asList(entry.getValue()));
		}
	}
	
	/**
	 * the array is interpreted as table in the following way for number of columns and rows: 
	 * <br>ar[numOfColumns][numOfRows], length of header String and number of columns must match
	 * <br>the columns are named Col0, Col1, Col2, ...
	 * @param ar
	 */
	public Table(Object[][] ar) {		
		for (int i = 0; i < ar.length; i++) {
			this.data.put("COL" + (i), Arrays.asList(ar[i]));
			++i;
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
		if (headers.length != ar.length) {
			throw new IllegalArgumentException("length of headers array (" + headers.length + ") and number of columns in ar (" + ar.length + ") must match");
		}
		int i = 0;
		for (String h : headers) {
			this.data.put(h, Arrays.asList(ar[i]));
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
	 * initializes the Table with empty Columns with headers specified
	 * @param headers
	 */
	public Table(String[] headers){
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
	 * get the column data as list of objects specified by columnName
	 * <br>this creates a deep clone of the underlying objects inside the list
	 * @param columnName
	 * @return
	 */
	public List<Object> getColumn(String columnName){
		List<Object> newList = new ArrayList<Object>();
		Collections.copy(this.data.get(columnName), newList);
		return newList;
	}
	
	/**
	 * return all columNames
	 * @return
	 */
	public String[] getColumnNames() {
		return (String[]) this.data.keySet().toArray();
	}
	
	/**
	 * return columnName at index c
	 * @param c
	 * @return
	 */
	public String getColumnName(int c) {
		this.checkRange(0, c);
		String[] columnNames = this.getColumnNames();
		return columnNames[c];
	}
	
	/**
	 * add empty column with columnName
	 * @param columnName
	 */
	public void addColumn(String columnName) {
		if (this.getNumberOfRows() == 0) {
			this.data.put(columnName, new ArrayList<Object>());
		} else {
			Object[] emptyColumn = new Object[this.getNumberOfRows()];
			this.addColumn(columnName, emptyColumn);
		}
	}

	/**
	 * add column with columnData
	 * @param columnData
	 */
	public void addColumn(Object[] columnData) {
		this.addColumn("COL" + (this.getNumberOfColumns()), columnData);
	}
	
	/**
	 * add a new column with columnName with columnData
	 * @param columnName
	 * @param columnData
	 */
	public void addColumn(String columnName, Object[] columnData) {
		if (this.hasColumn(columnName)) {
			throw new IllegalArgumentException("a column with name '" + columnName + "' is already present. Try replace or rename the column you want to add.");
		}
		this.data.put(columnName, Arrays.asList(columnData));
		this.fillUpColumns();
	}
	
	public void addColumn(List<Object> columnData) {
		String columnName = "COL" + this.getNumberOfColumns();
		this.addColumn(columnName, columnData);
	}
	
	public void addColumn(String columnName, List<Object> columnData) {
		if (this.hasColumn(columnName)) {
			throw new IllegalArgumentException("a column with name '" + columnName + "' is already present. Try replace or rename the column you want to add.");
		}
		this.data.put(columnName, columnData);
		this.fillUpColumns();
	}
	
	/**
	 * adds a new row based on the array rowData, if rowData does not equal available columns an Exception is thrown 
	 * @param rowData
	 */
	public void addRow(Object[] rowData) {
		if (rowData.length != this.getNumberOfColumns()) {
			throw new IllegalArgumentException("elements in rowData (" + rowData.length + ") does not equal the number of columns (" + this.getNumberOfColumns() + ") in this Table");
		}
		int c = 0;
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			entry.getValue().add(rowData[c]);
			++c;
		}
	}
	
	/**
	 * adds a new row based on the data in rowData, if rowData does not equal available columns an Exception is thrown
	 * @param rowData
	 */
	public void addRow(List<Object> rowData) {
		if (rowData.size() != this.getNumberOfColumns()) {
			throw new IllegalArgumentException("elements in rowData (" + rowData.size() + ") does not equal the number of columns (" + this.getNumberOfColumns() + ") in this Table");
		}
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
		this.data.get(columnName).add(value);
		this.fillUpColumns();
	}
	
	/**
	 * adds all columns from table to this, if overwrite is set to true,
	 * <br>then the column with the same name in table is kept
	 * @param table
	 */
	public void add(Table table, boolean overwrite) {
		for (Entry<String, List<Object>> entry : table.data.entrySet()) {
			this.data.put(entry.getKey(), entry.getValue());
		} 
	}
	
	/**
	 * adds all rows of table to the bottom of this object, if the number of columns is equal
	 * @param table
	 */
	public void append(Table table) {
		this.append(table, true);
	}
		
	/**
	 * Adds all rows of table to the bottom of this object, if the number of columns is equal.
	 * <br>If ignoreColumnNames is true, then only rows of columns from table present in this object are added,
	 * <br>else the rows are added to the same index of columns they came from in table
	 * @param table
	 * @param ignoreColumnNames
	 */
	public void append(Table table, boolean ignoreColumnNames) {
		// check if column size of table matches this table
		if (this.getNumberOfColumns() != table.getNumberOfColumns()) {
			throw new IllegalArgumentException("Number of columns do not match (" + this.getNumberOfColumns()  + " != " + table.getNumberOfColumns() + ")");
		}
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
	 * replace an element at column specified by columnName and row index with value
	 * @param columnName
	 * @param index
	 * @param value
	 */
	public void set(String columnName, int index, Object value) {
		this.checkRange(index, 0, columnName);
		this.data.get(columnName).set(index, value);
	}
	
	/**
	 * replace the whole row with new data
	 * @param row
	 * @param data
	 */
	public void set(int row, Object[] data) {
		if (data.length != this.getNumberOfColumns()) {
			throw new IllegalArgumentException("elements in rowData (" + data.length + ") does not equal the number of columns (" + this.getNumberOfColumns() + ") in this Table");
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
		this.checkRange(row, column);
		int c = 0;
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			if (c == column) {
				entry.getValue().set(row, value);
				return;
			}
		}
	}
	
	/**
	 * removes the value specified with columnName and row index
	 * @param columnName
	 * @param index
	 */
	public void remove(String columnName, int index) {
		// check if indices are within limits
		this.checkRange(index, 0, columnName);
		this.nullifyElement(columnName, index);
		if (this.isRowRemovable(index)) {
			this.removeRow(index);
		}
	}
	
	/**
	 * removes the whole row at index
	 * @param index
	 */
	public void remove(int index) {
		// check if indices are within limits
		this.checkRange(index);
		this.removeRow(index);
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
	 * removes a column specified by columnName
	 * @param columnName
	 */
	public void remove(String columnName) {
		this.data.remove(columnName);
	}
	
	/**
	 * removes the column at index
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
	 * removes duplicates in this Table based on every single column in otherTable
	 * @param otherTable
	 */
	public void removeDuplicates(Table otherTable) {
		this.removeDuplicates(otherTable, new int[0]);
	}
	
	public void removeDuplicates(Table otherTable, int[] columnIndices) {
		if (columnIndices.length != 0) {
			
		} else {
			
		}
	}
	
	public void removeDuplicates(Table otherTable, String[] columnNames) {
		
	}
	
	/**
	 * extracts all the columns specified by columnNames into new Table
	 * @param columnNames
	 * @return
	 */
	public Table extract(String[] columnNames) {
		Table t = new Table();	
		for (String columnName : columnNames) {
			// check if columnName is present in this table
			this.checkRange(0, 0, columnName);
			t.addColumn(columnName, this.getColumn(columnName));
		}	
		return t;
	}
	
	/**
	 * extracts all columns in this table specified by indices in columnIndices
	 * @param columns
	 * @return
	 */
	public Table extract(int[] columnIndices) {
		// check if all column indices are in range
		Table t = new Table();	
		for (int c : columnIndices) {
			this.checkRange(0, c);
			String columnName = this.getColumnName(c);
			t.addColumn(columnName, this.getColumn(columnName));
		}		
		return t;
	}
	
	public void sort(String sortColumn, SortMode sortMode) {
		
	}
	
	public Object xLookup(Object searchValue, String searchColumn, String lookupColumn, MatchMode matchMode, SearchMode searchMode) {
		Objects.requireNonNull(searchValue);
		this.checkColumn(searchColumn);
		this.checkColumn(lookupColumn);
		
		List<Object> searchData = new ArrayList<Object>();
		Collections.copy(this.data.get(searchColumn), searchData);
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
	
	public Table filter(Object filterValue, String filterColumn, Comparator comparator) {
		
		return null;
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
	public int getNumberOfColumns() {
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
		s[0] = this.getNumberOfRows();
		s[1] = this.getNumberOfColumns();
		return s;
	}
	
	/**
	 * returns the total number of elements in the table including NULL entries
	 * numOfElems = rows * cols
	 * @return
	 */
	public int getNumberOfElements() {
		return this.getNumberOfRows() * this.getNumberOfColumns();
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
	 * return the number of rows in Table
	 * @return
	 */
	public int getNumberOfRows() {
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			return entry.getValue().size();
		}
		return 0;
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
		if (this.getNumberOfRows() >= row) {
			return true;
		} else {
			return false;
		}
	}
	
	@Override
	public Table clone() {
		Table t = new Table();
		for (Entry<String, List<Object>> entry : this.data.entrySet()) {
			t.addColumn(entry.getKey(), entry.getValue());
		}
		return t;
	}
	
	/**
	 * checks whether column is within range in Table
	 * @param column
	 */
	public boolean isColumnIndexWithinLimit(int column) {
		if (this.getNumberOfColumns() >= column) {
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
			throw new IllegalArgumentException("column index (" + col + ") is out of range!");
		}
		if (!this.isRowIndexWithinLimit(row)) {
			throw new IllegalArgumentException("row index (" + row + ") is out of range!");
		}
		if (columnName != null) {
			if (!this.hasColumn(columnName)) {
				throw new IllegalArgumentException("Column (" + columnName + ") is not present in Table!");
			}
		}
	}
	
	private void checkColumn(String columnName) {
		if (!this.data.containsKey(columnName)) {
			throw new IllegalArgumentException(this.getClass().getSimpleName() + " does not contain column '" + columnName + "'.");
		}
	}

		
	/**
	 * string representation of this Table
	 * @return
	 */
	@Override
	public String toString() {
		String s = Table.class.getSimpleName() + "(" + this.getNumberOfRows() + ", " + this.getNumberOfColumns() + "):\n";
		String temp = this.data.toString();
		temp = temp.replace("], ", "]\n\t");
		temp = temp.replace("{", "{\n\t");
		temp = temp.replace("}", "\n}");
		return s + temp;
	}
	
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
	
}

