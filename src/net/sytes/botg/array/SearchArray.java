package net.sytes.botg.array;

import java.util.Arrays;
import java.util.List;

public class SearchArray {

	/**
	 * method returns true/false whether String str is contained in Array ar
	 * @param ar
	 * @param str
	 * @return true/false
	 */
	public static boolean isStringInArray(String[] ar, String str) {
		List<String> list = Arrays.asList(ar);
		if (list.contains(str)) {
			return true;
		} else {
			return false;
		}
	}			
	
    public static int countStringOccurencesInArray(String[] ar, String s, boolean caseSensitive) {
    	int count = -1;
        for (int i = 0; i < ar.length; i++) {
        	if (!caseSensitive) {
	        	if (s.toLowerCase().equals(ar[i].toLowerCase())){
	        		count++;
	            }
        	} else {
        		if (s.equals(ar[i])){
	        		count++;
	            }
        	}
        }
        return count;
    }  
    
    public static boolean isIntInArray(int[] ar, int i) {
    	for (int j = 0; j < ar.length; j++) {
    		if (ar[j] == i) {
    			return true;
    		}
    	}
    	return false;
    }
}
