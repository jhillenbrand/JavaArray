package datatypes.buffers;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import net.sytes.botg.datatypes.DataType;
import net.sytes.botg.datatypes.buffers.Buffer;
import net.sytes.botg.datatypes.buffers.TimedBuffer;

public class UnitTest_Buffer {

	@Test
	public void test000() {
		
		Buffer buf = new Buffer.Builder().build();
		
		System.out.println(buf);
		
	}
	
	@Test
	public void test010() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(10)
				.build();
		
		for (int i = 1; i < 15; i++) {
			buf.push(i);
			System.out.println(buf);
			System.out.println("newest:" + buf.newest());
			System.out.println("oldest:" + buf.oldest());
			System.out.println(Arrays.toString(buf.toArray(true)));
		}
		
		
	}
		
	@Test
	public void test020() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(3)
				.build();
				
		Object[] elements = new Object[]{1, 2, 3};
		
		buf.push(elements);
		
		System.out.println(buf.toString());
		System.out.println(buf.oldest());
		System.out.println(buf.newest());
	}
	
	@Test
	public void test021() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(3)
				.build();
				
		Object[] elements = new Object[]{1, 2};
		
		buf.push(elements);
		
		System.out.println(buf.toString());
		System.out.println(buf.oldest());
		System.out.println(buf.newest());
		
		buf.push(3);
		System.out.println(buf.toString());
	}
	
	@Test
	public void test030() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(3)
				.build();
				
		Object[] elements = new Object[]{2, 3};
		
		buf.push(1);
		
		buf.push(elements);
		
		System.out.println(buf.toString());		
		System.out.println(buf.oldest());
		System.out.println(buf.newest());
	}
	
	@Test
	public void test031() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(3)
				.build();
				
		Object[] elements = new Object[]{2, 3};
		
		buf.push(1);
		
		buf.push(elements);
		
		System.out.println(buf.toString());		
		System.out.println(buf.oldest());
		System.out.println(buf.newest());
		
		buf.push(4);
		System.out.println(buf.toString());
	}
	
	@Test
	public void test040() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(3)
				.build();
				
		Object[] elements = new Object[]{3, 4};
		
		buf.push(1);
		buf.push(2);
		
		buf.push(elements);
		
		System.out.println(buf.toString());
		System.out.println(buf.oldest());
		System.out.println(buf.newest());
	}
	
	@Test
	public void test041() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(3)
				.build();
				
		Object[] elements = new Object[]{3, 4};
		
		buf.push(1);
		buf.push(2);
		
		buf.push(elements);
		
		System.out.println(buf.toString());
		System.out.println(buf.oldest());
		System.out.println(buf.newest());
		
		buf.push(5);
		System.out.println(buf.toString());
	}
	
	@Test
	public void test050() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(3)
				.build();
				
		Object[] elements = new Object[]{1, 2, 3, 4};
		
		//buf.push(0);
		
		buf.push(elements);
		
		System.out.println(buf.toString());
		System.out.println(buf.oldest());
		System.out.println(buf.newest());
	}
	
	@Test
	public void test051() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(3)
				.build();
				
		Object[] elements = new Object[]{1, 2, 3, 4};
		
		//buf.push(0);
		
		buf.push(elements);
		
		System.out.println(buf.toString());
		System.out.println(buf.oldest());
		System.out.println(buf.newest());
		
		buf.push(5);
		System.out.println(buf.toString());
		
	}
	
	@Test
	public void test052() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(3)
				.build();
				
		Object[] elements = new Object[]{5, 6, 7, 8};
		
		//buf.push(0);
		
		buf.push(1);
		buf.push(2);
		buf.push(3);
		buf.push(4);
		
		buf.push(elements);
		
		System.out.println(buf.toString());
		System.out.println(buf.oldest());
		System.out.println(buf.newest());
		
		buf.push(9);
		System.out.println(buf.toString());
		System.out.println(buf.oldest());
		System.out.println(buf.newest());
	}
	
	/* NOTE comparison of push(Object[]) versions
	@Test
	public void test053() {
		
		int n = 1_000_000;
		
		Buffer buf = new Buffer.Builder()
				.capacity(5)
				.build();
		
		Object[] elements = new Object[]{5, 6, 7, 8, 9};
		
		long st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			buf.push(elements);
		}
		
		long et = System.nanoTime();
		
		double el = et - st;
		double sp = el / n; 
		
		System.out.println("Buffer push -> Sampling Period per Element [ns]: " + sp);
		
		buf.clear();
		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			buf.push2(elements);
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("Buffer push2 -> Sampling Period per Element [ns]: " + sp);
		
		buf.clear();
		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			buf.push(elements);
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("Buffer push -> Sampling Period per Element [ns]: " + sp);
		
		buf.clear();
		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			buf.push2(elements);
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("Buffer push2 -> Sampling Period per Element [ns]: " + sp);
		
	}	
	*/
	
	@Test
	public void test060() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(5)
				.build();
		
		Object[] elements = new Object[]{5, 6, 7, 8, 9};
		buf.push(elements);
		
		Object[] deqElems = buf.toArray(false, 3);
		
		System.out.println(buf.toString());
		System.out.println(Arrays.toString(deqElems));
		
		buf.push(10);
		System.out.println(buf.toString());
		
		System.out.println(buf.pull());
		System.out.println(buf.toString());
		
		
	}
	
	@Test
	public void test061() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(5)
				.build();
		
		Object[] elements = new Object[]{5, 6, 7, 8, 9};
		buf.push(elements);

		System.out.println(buf.toString(true));
		
		buf.push(10);
		buf.push(11);
		buf.push(12);
				
		System.out.println(buf.toString(true));
		
		Object[] deqElems = buf.toArray(false, 4);
		
		System.out.println(buf.toString(true));
		System.out.println(Arrays.toString(deqElems));
		
		buf.push(13);
		System.out.println(buf.toString(true));
		
		System.out.println(buf.pull());
		System.out.println(buf.toString(true));
		
	}
	
	@Test
	public void test062() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(5)
				.build();
		
		Object[] elements = new Object[]{5, 6, 7, 8, 9};
		buf.push(elements);
		
		Object[] ar = buf.toArray(false);
		
		System.out.println(buf.toString());
		
		System.out.println(Arrays.toString(ar));
		
	}
	
	@Test
	public void test070() {
		
		Buffer buf = new Buffer.Builder()
				.capacity(6)
				.dataType(DataType.INT)
				.build();
		
		int s = 8;
		
		int w = 0;
		int c = 0;
		
		while (true) {
			
			Integer[] samples = new Integer[s];
			
			for (int i = 0; i < s; i++) {
				samples[i] = c;
				++c;
			}
			
			buf.push(samples);
			
			++w;
			if (w > 15) {
				break;
			}
		}
		
		System.out.println(buf.toString());
		
	}
	
	@Test
	public void test080() {
		
		int c = 100;
		
		Buffer buf = new Buffer.Builder()
				.capacity(c)
				.dataType(DataType.INT)
				.build();
		
		for (int i = 0; i < c; i++) {
			
			buf.push(i);
			
		}
		
		System.out.println(buf.toString());
		
		buf.clear(10);
		
		System.out.println(buf.toString());
		
		double[] ar = new double[] {1.0, 2.0, 3.0};
		Object[] objs = new Object[ar.length];
		
		System.arraycopy(ar, 0, objs, 0, ar.length);
		
		System.out.println(Arrays.toString(objs));
		
	}
	
	@Test
	public void test090() throws InterruptedException {
		
		TimedBuffer tb = new TimedBuffer.Builder()
				.capacity(20)
				.id("B")
				.build();
		
		
		Thread t1 = new Thread(new PushTask(tb));
		Thread t2 = new Thread(new PullTask(tb));
		
		t1.start();
		t2.start();
		
		Thread.sleep(100_000);
		
	}	
	
	@Test
	public void test100() {
		TimedBuffer tb = new TimedBuffer.Builder()
				.capacity(3)
				.build();
		
		tb.push(0);
		tb.push(1);
		tb.push(2);
		tb.push(3);
		
		System.out.println(tb.pull().value());
		
	}
	
	private class PushTask implements Runnable {
		
		private TimedBuffer buffer = null;
		
		private int n = 0;
		
		public PushTask(TimedBuffer buffer){
			this.buffer = buffer;
		}
		
		@Override
		public void run() {
			
			while (true) {
				
				//for (int i = 0; i < 5; i++) {
					this.buffer.push(this.n);
					System.out.println("PushTask: " + Arrays.toString(this.buffer.toArray(true)));
					++this.n;
				//}
				
				try {
					Thread.sleep(20);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}				
			}
			
		}
		
	}
	
	private  class PullTask implements Runnable {
		
		private TimedBuffer buffer = null;
		
		public PullTask(TimedBuffer buffer){
			this.buffer = buffer;
		}
		
		@Override
		public void run() {
			while (true) {
				if (this.buffer.size() >= 5) {
					System.out.println("PullTask: " + Arrays.toString(this.buffer.toArray(false, 5)));
					try {
						Thread.sleep(100);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
	}
	
}
