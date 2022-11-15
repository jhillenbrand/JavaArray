package spectrum;
import java.util.Arrays;

import org.apache.commons.math3.complex.Complex;
import org.junit.Before;
import org.junit.Test;

import net.sytes.botg.featureextraction.FFT;
import net.sytes.botg.signals.sampled.SampledSine;

public class UnitTest_FFT {

	private double[] ar;
	
	@Before
	public void testData() {
		SampledSine sine = new SampledSine.Builder().sampleRate(100.0).a(1.0).f(5.0).n(0.2).p(0.0).build();		
		ar = sine.sample(1048576);		
	}
	
	@Test
	public void testFFT() {
		SampledSine sine = new SampledSine.Builder().sampleRate(100.0).a(1.0).f(5.0).n(0.1).p(0.0).build();		
		double[] ar = sine.sample(1024);
		
		Complex[] c = FFT.fft(ar, true); 
		
		System.out.println(Arrays.toString(c));
	}
	
	@Test
	public void testFFT2() {
		SampledSine sine = new SampledSine.Builder().sampleRate(100.0).a(1.0).f(5.0).n(0.1).p(0.0).build();		
		double[] ar = sine.sample(1024);
		
		double[] p = FFT.singleSidedSpectrum(ar); 
		
		System.out.println(Arrays.toString(p));
	}
	
	@Test
	public void testFFT3() {
		
		//Complex[] c = FFT.fft(ar);
		
		double[][] p = FFT.singleSidedSpectrum(ar, 100.0); 
		//Plot.line(ar);
		//Plot.line(p[0], p[1]);
		
		//System.out.println(Arrays.toString(p));
	}
	
	@Test
	public void testToComplex() {
		
		Complex[] c = FFT.toComplex(ar); 
		
	}
	
}
