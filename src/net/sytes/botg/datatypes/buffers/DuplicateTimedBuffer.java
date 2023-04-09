package net.sytes.botg.datatypes.buffers;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.UUID;

import net.sytes.botg.datatypes.ConfigOption;
import net.sytes.botg.datatypes.DataType;

public class DuplicateTimedBuffer extends TimedBuffer {
	
	@ConfigOption
	private List<String> duplicateIds;
	
	protected Map<String, TimedBuffer> duplicates;
	
	public DuplicateTimedBuffer() {
		this(new Builder());
	}
	
	private DuplicateTimedBuffer(Builder builder) {
		super(builder.capacity);
		this.id = builder.id;
		this.dataType = builder.dataType;
		this.description = builder.description;
		this.unit = builder.unit;
		this.duplicates = builder.duplicates;
		this.duplicateIds = builder.duplicateIds;
		if (this.duplicates != null) {
			List<String> previousDuplicateIds = new ArrayList<String>(this.duplicateIds);
			this.duplicateIds = new ArrayList<String>();
			for (Entry<String, TimedBuffer> entry : this.duplicates.entrySet()) {
				this.duplicateIds.add(entry.getKey());
			}
			if (!this.duplicateIds.isEmpty()) {
				if (!previousDuplicateIds.equals(this.duplicateIds)) {
					throw new IllegalArgumentException("Duplicates and duplicateIds are setup differently!");
				}
			}
		}
		if (this.duplicateIds != null && this.duplicates == null) {
			this.duplicates = new LinkedHashMap<String, TimedBuffer>();
			for (String duplicateId : this.duplicateIds) {
				TimedBuffer buffer = new TimedBuffer.Builder()
						.id(duplicateId)
						.capacity(this.capacity())
						.dataType(this.dataType)
						.description(this.description)
						.unit(this.unit)
						.build();
				
				this.duplicates.put(duplicateId, buffer);
			}
		} else {
			throw new IllegalArgumentException("Either duplicates or duplicateIds must be set!");
		}
	}

	public static class Builder {
		
		private String id = DuplicateTimedBuffer.class.getSimpleName() + " [" + UUID.randomUUID().toString() + "]";
		private DataType dataType = DataType.DOUBLE;
		private String description = null;
		private String unit = null;
		private int capacity = 1;
		private List<String> duplicateIds = null;
		private Map<String, TimedBuffer> duplicates = null;

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
		
		public Builder duplicateIds(List<String> duplicateIds) {
			this.duplicateIds = duplicateIds;
			return this;
		}
		
		public Builder duplicateId(String duplicateId) {
			this.duplicateIds = new ArrayList<String>();
			this.duplicateIds.add(duplicateId);
			return this;
		}
		
		public Builder duplicates(Map<String, TimedBuffer> duplicates) {
			this.duplicates = duplicates;
			return this;
		}

		public DuplicateTimedBuffer build() {
			return new DuplicateTimedBuffer(this);
		}
	}
		
	@Override
	public synchronized void push(long time, Object element) {
		 if (this.count == this.elems.length) {
         	this.dequeue();
         }
		 for (Entry<String, TimedBuffer> entry : this.duplicates.entrySet()) {
			 if (entry.getValue().count == entry.getValue().elems.length) {
				 entry.getValue().dequeue();
			 }
		 }
         this.enqueue(time, element);
         for (Entry<String, TimedBuffer> entry : this.duplicates.entrySet()) {
			 entry.getValue().enqueue(time, element);
		 }
	}
	
	@Override
	public synchronized void push(long[] times, Object[] elements) {
		this.enqueue(times, elements);
		for (Entry<String, TimedBuffer> entry : this.duplicates.entrySet()) {
			 entry.getValue().enqueue(times, elements);
		}
	}

	@Override
	public synchronized void push(Object element) {
		 if (this.count == this.elems.length) {
         	this.dequeue();
         }
		 for (Entry<String, TimedBuffer> entry : this.duplicates.entrySet()) {
			 if (entry.getValue().count == entry.getValue().elems.length) {
				 entry.getValue().dequeue();
			 }
		 }
         this.enqueue(element);
         for (Entry<String, TimedBuffer> entry : this.duplicates.entrySet()) {
			 entry.getValue().enqueue(element);
		 }
	}
	
	@Override
	public synchronized void push(Object[] elements) {
		this.enqueue(elements);
		for (Entry<String, TimedBuffer> entry : this.duplicates.entrySet()) {
			 entry.getValue().enqueue(elements);
		 }
	}
	
	/**
    * returns a String describing content and config
    */
	@Override
	public String toString() {		
		StringBuilder sb = new StringBuilder();
		sb.append(this.getClass().getSimpleName()).append(" id=[").append(this.id).append("]");    	
		sb.append("\n\tElements:" + this.size() + "/" + this.capacity());
		sb.append("\n\tDataType: " + this.dataType);
		sb.append("\n\tDuplicateIds: " + this.duplicateIds);
		sb.append("\n\t" + unsynchronizedToTimeSeries());
		return sb.toString();    	
	}
	
	public List<String> getDuplicateIds() {
		return this.duplicateIds;
	}

	public Map<String, TimedBuffer> getDuplicates() {
		return this.duplicates;
	}
}
