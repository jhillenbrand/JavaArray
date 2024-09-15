package spectrum;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.math.Scalar;
import net.sytes.botg.array.math.Vec;
import net.sytes.botg.array.spectrum.Complex;
import net.sytes.botg.array.spectrum.Spectrum;

public class UnitTest_Spectrum {

	@Test
	public void test000() {
		
		int n = 10;
		
		double[] t = Vec.linspace(0, 4 * 2 * Math.PI, 2048);
		double[] s = new double[t.length];
		for (int i = 0; i < s.length; i++) {
			s[i] = Math.cos(5 * t[i]);
		}
		
		long st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			
			Spectrum.dft(s, true);
			
		}
		
		long et = System.nanoTime();
		
		double el = et - st;
		double sp = el / n; 
		
		System.out.println("DFT (pos) -> Sampling Period per Iteration [ns]: " + sp);
		
		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			
			Spectrum.dft(s, false);
			
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("DFT -> Sampling Period per Iteration [ns]: " + sp);
		
		st = System.nanoTime();
				
		for (int i = 0; i < n; i++) {
			
			Spectrum.singleSidedSpectrum(s);
			
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("FFT -> Sampling Period per Iteration [ns]: " + sp);
	}
	
	@Test
	public void test010() {
		
		double[] d = Vec.rand(10);
		
		Complex[] c = Spectrum.fft2(d);
		
		System.out.println(Arrays.toString(c));
		
	}
	
	@Test
	public void test020() {
		
		double f = 50.0;
		double fs = 1000.0;
		
		double[] sine = Vec.sine(1024, 1.0, f, fs);
		
		double[][] P = Spectrum.singleSidedSpectrum(sine, fs);
		
		double mf = Spectrum.meanFrequency(sine, fs);
		
		assertEquals(Scalar.round(mf), (int) f);
		
	}
	
}
