package net.sytes.botg.datatypes;

public class Unsigned {

	public static long uint(int i) {
		long li = i & 0xffffffffL;
		return li;
	}
	
	public static int sint(long lo) {
		Long lon = (Long) lo;
		return lon.intValue();
	}
	
	public static int ushort(short s) {
		int i = s & 0xffff;
		return i;
	}
	
	public static short sshort(int i) {
		Integer in = (Integer) i;
		return in.shortValue();
	}
	
	public static int ubyte(byte b) {
		int i = b & 0xff;
		return i;		
	}
	
	public static int[] ubyte(byte[] bs) {
		int[] ints = new int[bs.length];
		for (int s = 0; s < ints.length; s++) {
			ints[s] = ubyte(bs[s]);
		}
		return ints;
	}
	
}
