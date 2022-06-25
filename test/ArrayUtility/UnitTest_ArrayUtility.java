package ArrayUtility;

import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;
import org.junit.Test;

import net.sytes.botg.array.ArrayUtility;
import net.sytes.botg.array.math.Vec2Vec;

public class UnitTest_ArrayUtility {

	private static final int SIZE = 10_000_000;
	
	@Test
	public void test00() {
		
		System.out.println("WARMUP");
		System.out.println(Integer.MAX_VALUE);
	}	
	
	@Test
	public void test01() {
		
		ArrayUtility.createRandomDoubleArray(SIZE);
		
	}
	
	@Test
	public void test02() {
		
		Double[] ar = {10.0, 128.234, 3984.123};
		Double d = 10.9877;
		
		System.out.println(Arrays.toString(ar));
		System.out.println(d);
		
		Double[] nAr = Vec2Vec.append(ar, d);
		System.out.println(Arrays.toString(nAr));
	}
	
	@Test
	public void test03() {
		
		double[] ar = {10.0, 128.234, 3984.123};
		double d = 10.9877;
		
		System.out.println(Arrays.toString(ar));
		System.out.println(d);
		
		double[] nAr = Vec2Vec.append(ar, d);
		System.out.println(Arrays.toString(nAr));
	}
	
	@Test
	public void test040() {
		
		Character[] cw = new Character[2];
		cw[0] = 'A';
		cw[1] = 'B';
		
		char[] c = ArrayUtils.toPrimitive(cw);
		
		System.out.println(c);		
	}
	
	@Test
	public void test050() {
		
		double[] d = ArrayUtility.createRandomDoubleArray(100);
		
		double[] dd = ArrayUtility.subArray(d, 0, 100);
		
		System.out.println(Arrays.toString(dd));
		
	}
	
	@Test
	public void test060() {
		
		double[] d = ArrayUtility.createRandomDoubleArray(100);
		
		double[] dd = ArrayUtility.subArray(d, 1, 90);
		
		System.out.println(Arrays.toString(dd));
		
	}

}
