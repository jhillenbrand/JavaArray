package net.sytes.botg.array.test;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.ConvertArray;

public class UnitTest_ArrayUtility {

	
	private static final String FOLDER = "C:\\Users\\aeadmin\\Downloads\\test\\";
	
	@Test
	public void testNumberToByte() {
		
		//Number n = 0;
		double n = 1020480.56485;
		
		byte[] bytes = ConvertArray.convertDoubleToBytes2(n, 64);
		//byte[] bytes = ArrayUtility.convertIntToBytes(n, 16);
		System.out.println(Arrays.toString(bytes));
		
		FileOutputStream fos;
		try {
			fos = new FileOutputStream(FOLDER + "test1.bin");
			
			fos.write(bytes);
			
			fos.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	@Test
	public void testNumberArrayToByteFile() {
		Number[] ns = {120.2, 2134123, -34, 12};
		
		FileOutputStream fos;
		try {
			fos = new FileOutputStream(FOLDER + "test3.bin");
			byte[] bytes = null;
			for (Number n : ns) {
				bytes = ConvertArray.convertNumberToBytes(n);
				fos.write(bytes);
			}
			
			
			fos.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	@Test
	public void testShortToByte() {
		
		Short n = 32767;
		//double n = 1020480.56485;
		
		byte[] bytes = ConvertArray.convertShortToBytes(n);
		//byte[] bytes = ArrayUtility.convertIntToBytes(n, 16);
		System.out.println(Arrays.toString(bytes));
		
		FileOutputStream fos;
		try {
			fos = new FileOutputStream(FOLDER + "test2.bin");
			
			fos.write(bytes);
			
			fos.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
}
