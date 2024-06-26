package math;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.Test;
import org.junit.jupiter.api.DisplayName;

import net.sytes.botg.array.Ar;
import net.sytes.botg.array.math.Vec;
import net.sytes.botg.array.math.Vec.GroupBy;

public class UnitTest_Vec {

	double[] ar = {100, Double.NaN, 129812.9, 120.2174914, Double.NaN};
	double[] ar1 = Vec.rand(1000);
	double[] ar2 = Vec.rand(998);
	
	@Test
	public void test000() {
		
		double[] ar = {1.8, 2.0, 8.1, 9.0, 243.2, 2134.123, 23.0, -1.0, 12.9};
		
		double[] newAr = Vec.windowedSpan(ar, 3, true);
		
		System.out.println(Arrays.toString(newAr));
		
	}
	
	@Test
	public void test010() {
		
		double[] ar = {1.8, 2.0, 8.1, 9.0, 243.2, 2134.123, 23.0, -1.0, 12.9, -5.0, 4.5};
		
		double[] newAr = Vec.windowedSpan(ar, 3, false);
		
		System.out.println(Arrays.toString(newAr));
		
	}
	
	@Test
	public void test020() {
		
		double[] x = {22.2,	22.25,	22.3,	22.25,	22.2,	22.3,	22.2,	22.4,	22.4};
		
		double[] xUp = Vec.upsample2(x, 3);
		
		System.out.println(Arrays.toString(xUp));
	}
	
	@Test
	public void test030() {
		
		double[] x = {22.2,	22.25,	22.3,	22.25,	22.2,	22.3,	22.2,	22.4,	22.4};
		
		double[] y = Vec.zeroRange(x, 2, 5);
		
		System.out.println(Arrays.toString(x));
		System.out.println(Arrays.toString(y));
	}
	
	@Test
	public void test040() {
		double[] ar = Vec.linspace(0.0, 10.0, 100);
		
		double[] ar2 = Vec.downsampleBruteForce(ar, 100);
		
		System.out.println(Arrays.toString(ar2));
	}
	
	@Test
	public void test041() {
		double[] ar = Vec.linspace(0.0, 10.0, 100);
		
		double[] ar2 = Vec.downsampleBruteForce(ar, 50);
		System.out.println(Arrays.toString(ar2));
	}
	
	@Test
	public void test042() {
		double[] ar = Vec.linspace(0.0, 10.0, 100);
		
		double[] ar2 = Vec.downsampleBruteForce(ar, 33);
		System.out.println(Arrays.toString(ar2));
	}
	
	@Test
	public void test043() {
		
		int f = 11;
		
		double[] ar = Vec.linspace(0.0, 10.0, 100);
		
		double[] ar2 = Vec.downsampleBruteForce(ar, ar.length / f);
		
		System.out.println(Arrays.toString(ar2));
	}
	
	@Test
	public void test044() {
		
		int m = 5;
		
		double[] ar = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
		
		double[] ar2 = Vec.downsampleMax(ar, m);
		
		System.out.println(Arrays.toString(ar2));
	}
		
	@Test
	public void test050() {
		double[] ar = {0.0, 1.0, 2.0, 1.0, 0.0, -1.0, -1.0, 1.0, 0.0, 1.0};
		
		int[] inds = Vec.zeroCrossings(ar);
		
		System.out.println(Arrays.toString(inds));
		
	}
	
	@Test
	public void test051() {
		double[] ar = {0.0, 1.0, 2.0, 1.0, 0.0, -1.0, -1.0, 1.0, 0.0, 1.0};
		
		int[] inds = Vec.zeroCrossings(ar, 3);
		
		System.out.println(Arrays.toString(inds));
		
	}

	@Test
	public void test052() {
		
		double[] d = Vec.rand(100);
		
		double[] dd = Ar.sub(d, 0, 50);
		
		System.out.println(Arrays.toString(dd));
		
	}
	
	@Test
	public void test053() {
		
		double[] d = Vec.rand(100);
		
		double[] dd = Ar.sub(d, 1, 90);
		
		System.out.println(Arrays.toString(dd));
		
	}
	
	@Test
	public void test060() {
		
		double[] ar = {0.0, 0.0, -1.0, -0.5, 0.0, 1.5, 2.0, 2.4, 2.3, 2.3, 2.8, 1.0};
		
		int[] extrema = Vec.findLocalExtrema(ar);
		
		System.out.println(Arrays.toString(extrema));
		
	}
	
	@Test
	public void test070() {
		//System.out.println(MathArray.sum2(ar2));
		Vec.sum2(ar2);
	}
	
	
	@Test
	public void test071() {
		//System.out.println(MathArray.sum(ar2));
		Vec.sum(ar2);
	}
	
	@Test
	public void test080() {
		double[] x_s = {0, 1, 2, 3, 4};
		double[] y_s = {21,24,24,18,16};
		
		double[] x = {0,0.5,1,1.5,2,2.5,3,3.5,4};
		
		double[] y = Vec.splineInterp(x_s, y_s, x);
		
		System.out.println(y);
		
	}
	
	@Test
	public void test090() {
		double[] nans = Vec.nan(10);
		
		System.out.println(Arrays.toString(nans));
		
	}	

	@Test
	public void test100() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
		int numSamples = x.length;
		int fftSampleSize = 4;
		int overlapFactor = 2;
		int a = 0;
		
		if (overlapFactor > 1) {

			int numOverlappedSamples = numSamples * overlapFactor;
			int backSamples = fftSampleSize * (overlapFactor - 1) / overlapFactor;
			int fftSampleSize_1 = fftSampleSize - 1;
			double[] overlapAmp = new double[numOverlappedSamples];
			a = 0;
			for (int i = 0; i < numSamples; i++) {
				overlapAmp[a++] = x[i];
				if (a % fftSampleSize == fftSampleSize_1) {
					// overlap
					i -= backSamples;
				}
			}
			numSamples = numOverlappedSamples;
			x = overlapAmp;
		}
		
		System.out.println(Arrays.toString(x));
		
	}
	
	@Test
	public void test110() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
		
		double[][] W = Vec.overlapWindows(x, 4, 2);
		
		Ar.print(W);
		
	}	
	
	@Test
	public void test111() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] W = Vec.overlapWindows(x, 4, 2);
		
		Ar.print(W);
		
	}	
	
	@Test
	public void test112() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] W = Vec.overlapWindows(x, 4, 1);
		
		Ar.print(W);
		
	}	
	
	@Test
	public void test113() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] W = Vec.overlapWindows(x, 3, 2);
		
		Ar.print(W);
		
	}	
	
	@Test
	public void test114() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] W = Vec.overlapWindows(x, 3, 1);
		
		Ar.print(W);
		
	}	
	
	@Test
	public void test115() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] W = Vec.overlapWindows(x, 5, 2);
		
		Ar.print(W);
		
	}	
	
	@Test
	public void test116() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] W = Vec.overlapWindows(x, 5, 3);
		
		Ar.print(W);
		
	}	
	
	@Test
	public void test117() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] W = Vec.overlapWindows(x, 5, 0);
		
		Ar.print(W);
		
	}	
	
	@Test
	public void test118() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
		
		double[][] W = Vec.overlapWindows(x, 5, 4);
		
		Ar.print(W);
		
	}	
	
	@Test
	public void test119() {
	
		double[] x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
		
		double[][] w = Vec.overlapWindows(x, 3, 2);
		
		Ar.print(w);
		
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 2")
	public void test121() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(Vec.isMonotone(ar, false), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 3")
	public void test122() {
		double[] ar = {1.0, 1.0, 3.0, 4.0, 4.1};
		
		assertEquals(Vec.isMonotone(ar, true), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 4")
	public void test123() {
		double[] ar = {1.0, 0.99, 3.0, 4.0, 4.1};
		
		assertEquals(Vec.isMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 5")
	public void test124() {
		double[] ar = {1.0, 0.99, 0.9, 0.0, -1.1};
		
		assertEquals(Vec.isMonotone(ar, false), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 1")
	public void test125() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(Vec.isStrictMonotone(ar, true), true);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 2")
	public void test126() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(Vec.isStrictMonotone(ar, false), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 3")
	public void test127() {
		double[] ar = {1.0, 1.0, 3.0, 4.0, 4.1};
		
		assertEquals(Vec.isStrictMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 4")
	public void test128() {
		double[] ar = {1.0, 0.99, 3.0, 4.0, 4.1};
		
		assertEquals(Vec.isStrictMonotone(ar, true), false);
	}
	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isStrictMonotone() - 5")
	public void test129() {
		double[] ar = {1.0, 0.99, 0.9, 0.0, -1.1};
		
		assertEquals(Vec.isStrictMonotone(ar, false), true);
	}
	
	@Test
	public void test130() {
		
		double[] x = {0, 1, 2, 3, 4};
		double[] y = {21, 24, 24, 18, 16};
		
		double[] x_s = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4};
		
		double[] y_s = Vec.splineInterp(x, y, x_s);
		
		Ar.print(y_s);
		
	}
	
	@Test
	public void test132() {
		
		double[] logspace = Vec.logspace(1, 10, 10);
		System.out.println(Arrays.toString(logspace));
		
	}
	

	
	@Test
	@DisplayName("Testing Monotonicity Methods -> isMonotone() - 1")
	public void test131() {
		double[] ar = {1.0, 2.0, 3.0, 4.0, 4.1};
		
		assertEquals(Vec.isMonotone(ar, true), true);
	}

	@Test
	public void test133() {
		double[] nans = Vec.nan(10);
		double[] rands = Vec.rand(10);
		
		double[] prods = Vec.product(rands, nans);
		
		System.out.println(Arrays.toString(prods));
		
	}
	
	@Test
	public void test140() {
		
		double[] x = new double[] {1.0, 2.0, 3.0};
		
		double[] y = Ar.copy(x);
		
		y[0] = 0.0;
		
		System.out.println(Arrays.toString(x));
		
	}

	
	@Test
	public void test151() {
		
		int s = 1_000;
		int n = 100_000;		

		//int[] ints = ArUtils.createRandomIntArray(s);
		int[] ints = Vec.linspace(1, 100);
				
		long t1 = System.nanoTime();
		
		
		for (int i = 0; i < ints.length; i++) {
			Ar.checkForGreaterEqualZero2(ints);
		}
		
		long t2 = System.nanoTime();
		long t_e = t2 - t1;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		
		
		
		t1 = System.nanoTime();		
		for (int i = 0; i < ints.length; i++) {
			Ar.checkForGreaterEqualZero(ints);
		}
		
		t2 = System.nanoTime();
		t_e = t2 - t1;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
	}
	
	@Test
	public void test160() {
		
		double[] values = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
		int[] index = {0, 1, 0, 2, 3, 1, 2, 0, 3};
		double[][] groups = Vec.group(values, index);
		
		Ar.print(groups);
		
	}

	@Test
	public void test201() {
		double[] data = Vec.rand(10000);
		
	}
	
	@Test
	public void test202() {
		
		Vec.rand(1000);
		
	}
	
	@Test
	public void test203() {
		
		Double[] ar = {10.0, 128.234, 3984.123};
		Double d = 10.9877;
		
		System.out.println(Arrays.toString(ar));
		System.out.println(d);
		
		Double[] nAr = Ar.append(ar, d);
		System.out.println(Arrays.toString(nAr));
	}
	
	@Test
	public void test204() {
		
		double[] ar = {10.0, 128.234, 3984.123};
		double d = 10.9877;
		
		System.out.println(Arrays.toString(ar));
		System.out.println(d);
		
		double[] nAr = Ar.append(ar, d);
		System.out.println(Arrays.toString(nAr));
	}
		
	@Test
	public void test301() {
		
		System.out.println(Arrays.toString(Vec.isNaN(ar)));
		
	}
	
	@Test
	public void test302() {
		
		double[][] dd = Vec.separateDataIntoWindows(ar1, 100, true);
		
		System.out.println(Arrays.toString(dd));
		
	}
	
	@Test
	public void test303() {
		
		double[][] dd = Vec.separateDataIntoWindows(ar2, 100, false);
		
		System.out.println(Arrays.toString(dd));
		
	}
	
	@Test
	public void test304() {
		
		double[] x = {1.0, 1.1, 9.1, 76.3, 1.1, 2.0, 9.1};
		
		double[] xu = Ar.unique(x);
		
		System.out.println(Arrays.toString(xu));
		
	}
	
	@Test
	public void test305() {
		
		double[] x = {1.0, 1.1, 9.1, 76.3, 1.1, 2.0, 9.1};
		
		double[] xu = Ar.unique(x);
		
		System.out.println(Arrays.toString(xu));
		
	}
	
	@Test
	public void test306() {
		
		double[] x = {1.0, 1.1, 9.1, 76.3, 1.1, 2.0, 9.1};
		
		int n = 1_000_000;
		long t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			double[] xu = Ar.unique(x);
			//int b = bs.length;
		}
		long t2 = System.nanoTime();
		long t_e = t2 - t1;
		
		System.out.println("unique2: " + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		
		
		t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			double[] xu = Ar.unique(x);
			//int b = bs.length;
		}
		t2 = System.nanoTime();
		t_e = t2 - t1;
		
		System.out.println("unique: " + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
				
	}
	
	@Test
	public void test400() {
		byte[] bs = {0, 1, 0, 1, 0, 1, 0, 1};
		
		int n = 10_000;
		long ns = 1_000_000_000;
		long t1 = System.nanoTime();
		for (int i = 0; i < n; i++) {
			Vec.oddToEven(bs);
			//int b = bs.length;
		}
		long t2 = System.nanoTime();
		long t_e = t2 - t1;
		double sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
		
		for (int i = 0; i < n; i++) {
			Vec.oddToEven2(bs);
			//int b = bs.length;
		}
		
		t2 = System.nanoTime();
		t_e = t2 - t1;
		sr = (double) ns / t_e * n;
		
		System.out.println("" + (t2 - t1) + " [ns]; " + (((double) (t2 - t1)) / (double) n) + "[ns]");
		System.out.println("sr=" + sr + "[Hz]");
	}
	
	@Test
	public void test401() {
		
		double[] ar1 = {1, 2, 3, 4, 5}; 
		System.out.println(Arrays.toString(Vec.even(ar1)));
		
		System.out.println(Arrays.toString(Vec.odd(ar1)));
		
	}
	
	@Test
	public void test410() {
		
		int[] r = Vec.randInt(10, 1, 10);
		
		int[] inds = Vec.bubbleSortInd(r, true);
		
		Ar.print(r);
		Ar.print(inds);
		
		Ar.print(Ar.elementsAt(r, inds));
		
	}
	
	@Test
	public void test420() {
		
		double[] x = Vec.linspace(5);
		double[] y = Vec.rand(5);
		
		double[] y_reg = Vec.linReg(x, y);
		
		Ar.print(y);
		Ar.print(y_reg);
		
	}
	
	@Test
	public void test430() {
		
		double[] x = {1, 1};
		double[] xn = Vec.unitVector(x);
		
		Ar.print(xn);
		
	}
	
	@Test
	public void test440() {
		
		double[] x = Vec.linspace(10);
		
		double[] y = Vec.padding(x, 15);
		
		Ar.print(x);
		Ar.print(y);
		
	}
	
	@Test
	public void test441() {
		
		double[] x = Vec.linspace(10);
		
		double[] y = Vec.padding(x, 15, 1.0);
		
		Ar.print(x);
		Ar.print(y);
		
	}
	
	@Test
	public void test442() {
		
		double[] x = Vec.linspace(5);
		
		double[] y = Vec.mirroredPadding(x, 8);
		
		Ar.print(x);
		Ar.print(y);
		
	}
	
	@Test
	public void test443() {
		
		double[] x = Vec.linspace(5);
		
		double[] y = Vec.mirroredPadding(x, 25);
		
		Ar.print(x);
		Ar.print(y);
		
	}
	
	
	@Test
	public void test450() {
		double[] x = Vec.linspace(9);
		double[][] X = Vec.matrix(x, 3, 3, false);
		
		Ar.print(X);
	}
	
	@Test
	public void test451() {
		
		double[] x = {1, 1, 1};
		
		double[][] X = Vec.matrix(x, 1);
		
		Ar.print(X);
		
	}
	
	@Test
	public void test452() {
		
		double[] x = {1, 1, 1};
		
		double[][] X = Vec.matrix(x, 3, 1, true);
		
		Ar.print(X);
		
	}
	
	@Test
	public void test460() {
		
		double[] x = Vec.rand(20, -5.0, 5.0);
		
		double[] y = Vec.scale(x, -1.0, 1.0);
		
		double[] z = Vec.scale(x, 0.0, 1.0);
		
		Ar.print(y);
		Ar.print(z);
		
	}
	
	@Test
	public void test470() {
		
		double[] x = Vec.rand(10, 0.0, 100.0);
		
		double[] max = Vec.maxk(x, 3);
		int[] i = Vec.maxkInd(x, 3);
		
		Ar.print(x);
		Ar.print(max);
		Ar.print(i);
		
	}
	
	@Test
	public void test471() {
		
		double[] x = Vec.rand(10, 0.0, 100.0);
		
		double[] min = Vec.mink(x, 3);
		int[] i = Vec.minkInd(x, 3);
		
		Ar.print(x);
		Ar.print(min);
		Ar.print(i);
		
	}
	
	@Test
	public void test480() {
		
		double[] x = Vec.linspace(100);
		double[] y = Vec.rand(100, 0.0, 100.0);
		
		double[][] d = Vec.downsampleNonUniformDataByDeletingCloseGradients(x, y, 80);
		
		Ar.print(d[0]);
		Ar.print(d[1]);
		
	}
	

	
	@Test
	public void test490() {
		
		double[] sourceValues = new double[]{0.0, 1.0, 2.0, 3.0};
		
		double[] searchValues = new double[]{0.5, 1.1, 1.9, 2.0};
		double[] returnValues = new double[]{2.5, 5.0, 10.0, 20.0};
		
		double[] mappedValues = Vec.mapValues(sourceValues, searchValues, returnValues);
		
		Ar.print(mappedValues);
		
	}
	
	@Test
	public void test500() {
		double[] x = Vec.linspace(10);
		
		double d = 3.0;
		
		Ar.print(x);
		
		System.out.println(Vec.findClosest(x, d));
		
	}
	
	@Test
	public void test501() {
		double[] x = Vec.rand(10, 0.0, 10.0);
		
		double d = 3.0;
		
		Ar.print(x);
		
		System.out.println(Vec.findClosest(x, d));
		
	}
	
	@Test
	public void test510(){
		
		double[] x = {1.0, 2.0, 1.0, 3.0, 1.0, 2.0, 2.0, 3.0};
		double[] y = {3.0, 4.0, 2.0, 10.0, 5.5, 3.0, 2.1, 9.0};
		
		double[][] z = Vec.group(x, y, GroupBy.MEAN);
		
		Ar.print(z);
		
	}
	
	@Test
	public void test520() {
		
		double[] x = Vec.linspace(3);
		
		Ar.print(x);
		
		double[][] X = Vec.transpose(x);
		
		Ar.print(X);
		
	}
	

	
	@Test
	public void test530() {
		double[] x = Vec.linspace(10);
		
		double[] y = Vec.linspace(5);
		
		double[] z = Vec.append(x, y);
		
		Ar.print(z);
	}
	
	@Test
	public void test531() {
		double[] a = Vec.linspace(10);
		
		double[] b = Vec.linspace(5);
		
		double[] c = Vec.linspace(3);
		
		double[] d = Vec.append(a, b, c);
		
		Ar.print(d);
	}
	
	@Test
	public void test532() {
		double[] a = {};
		
		double[] b = Vec.linspace(10);
		
		double[] c = Vec.linspace(3);
		
		double[] d = Vec.append(a, b, c);
		
		Ar.print(d);
	}
		
}
