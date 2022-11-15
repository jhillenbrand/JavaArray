package spectrum;

import org.junit.jupiter.api.Test;

import net.sytes.botg.featureextraction.DFT;
import net.sytes.botg.featureextraction.FFT;
import net.sytes.botg.plotting.SubPlot;
import net.sytes.botg.plotutils.P;
import net.sytes.botg.signals.sampled.SampledSine;

public class UnitTest_DFT {

	@Test
	public void test000() {
		
		int n = 10;
		
		SampledSine sine = new SampledSine();
		double[] data = sine.sample(2048);		
		
		long st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			
			DFT.dft2(data, true);
			
		}
		
		long et = System.nanoTime();
		
		double el = et - st;
		double sp = el / n; 
		
		System.out.println("DFT2 (pos) -> Sampling Period per Iteration [ns]: " + sp);
		
		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			
			DFT.dft2(data, false);
			
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("DFT2 -> Sampling Period per Iteration [ns]: " + sp);
		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			
			DFT.dft(data, false);
			
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("DFT (pos) -> Sampling Period per Iteration [ns]: " + sp);

		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			
			DFT.dft(data, false);
			
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("DFT -> Sampling Period per Iteration [ns]: " + sp);

		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			
			DFT.dft3(data, true);
			
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("DFT3 (pos) -> Sampling Period per Iteration [ns]: " + sp);
		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			
			DFT.dft3(data, false);
			
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("DFT3 -> Sampling Period per Iteration [ns]: " + sp);
		
		st = System.nanoTime();
		
		for (int i = 0; i < n; i++) {
			
			FFT.singleSidedSpectrum(data);
			
		}
		
		et = System.nanoTime();
		
		el = et - st;
		sp = el / n; 
		
		System.out.println("FFT -> Sampling Period per Iteration [ns]: " + sp);
	}
	
	@Test
	public void test010() {
		
		SampledSine sine = new SampledSine.Builder()
				.sampleRate(100.0)
				.build();
		
		double[] data = sine.sample(2048);	
	
		SubPlot sp = new SubPlot(7);
		
		double[] p = null;
		
		sp.sub(0);
		sp.title("FFT");
		p = FFT.singleSidedSpectrum(data);
		sp.addLine(p);
		sp.drawnow();
		
		sp.sub(1);
		sp.title("DFT");
		p = DFT.dft(data, true);
		sp.addLine(p);
		sp.drawnow();
		
		sp.sub(2);
		sp.title("DFT2");
		p = DFT.dft2(data, true);
		sp.addLine(p);
		sp.drawnow();
		
		sp.sub(3);
		sp.title("DFT3");
		p = DFT.dft3(data, true);
		sp.addLine(p);
		sp.drawnow();
		
		sp.sub(4);
		sp.title("DFT3_");
		p = DFT.dft3(data, false);
		sp.addLine(p);
		sp.drawnow();
		
		sp.sub(5);
		sp.title("DFT2_");
		p = DFT.dft2(data, false);
		sp.addLine(p);
		sp.drawnow();
		
		sp.sub(6);
		sp.title("DFT_");
		p = DFT.dft(data, false);
		sp.addLine(p);
		sp.drawnow();
		
		P.pause();
	}
	
}
