package spectrum;

import org.junit.jupiter.api.Test;

import net.sytes.botg.array.math.Vec;
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
	
}
