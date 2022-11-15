package net.sytes.botg.array.spectrum;

import org.apache.commons.math3.util.FastMath;

import com.github.psambit9791.jdsp.transform.DiscreteFourier;

public class DFT {
			
	/**
	 * returns the absolute values of the amplitude spectrum of {@code x},
	 * <br>if {@code onlyPositive} is set to true, the only the non-mirrored output is returned
	 * @param x
	 * @param onlyPositive
	 * @return
	 */
	public static double[] dft(double[] x, boolean onlyPositive) {
		int n = x.length;
		double[] p;
		if (onlyPositive) {
			p = new double[n / 2];
	        for (int i = 0; i < n / 2; i++) {
	            double real = 0;
	            double imag = 0;
	            for (int j = 0; j < n; j++) {
	                double angle = (2 * Math.PI * i * j) / n;
	                real += x[j] * Math.cos(angle);
	                imag += x[j] * Math.sin(angle);
	            }
	            p[i] = Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
	        }
		} else {
			p = new double[n];
	        for (int i = 0; i < n; i++) {
	            double real = 0;
	            double imag = 0;
	            for (int j = 0; j < n; j++) {
	                double angle = (2 * Math.PI * i * j) / n;
	                real += x[j] * Math.cos(angle);
	                imag += x[j] * Math.sin(angle);
	            }
	            p[i] = Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
	        }
		}
        return p;
	}

	/**
	 * returns the absolute values of the amplitude spectrum of {@code x},
	 * <br>if {@code onlyPositive} is set to true, the only the non-mirrored output is returned
	 * <br>uses com.github.psambit9791.jdsp.transform.DiscreteFourier DFT implementation
	 * @param x
	 * @param onlyPositive
	 * @return
	 */
	public static double[] dft2(double[] x, boolean onlyPositive) {
		DiscreteFourier fft = new DiscreteFourier(x);
		fft.dft();
		double[] p = fft.returnAbsolute(onlyPositive);
		return p;
	}
	
	/**
	 * returns the absolute values of the amplitude spectrum of {@code x},
	 * <br>if {@code onlyPositive} is set to true, the only the non-mirrored output is returned
	 * <br>uses the apache commons math fast api
	 * @param x
	 * @param onlyPositive
	 * @return
	 */
	public static double[] dft3(double[] x, boolean onlyPositive) {
		int n = x.length;
		double[] p;
		if (onlyPositive) {
			p = new double[n / 2];
	        for (int i = 0; i < n / 2; i++) {
	            double real = 0;
	            double imag = 0;
	            for (int j = 0; j < n; j++) {
	                double angle = (2 * Math.PI * i * j) / n;
	                real += x[j] * FastMath.cos(angle);
	                imag += x[j] * FastMath.sin(angle);
	            }
	            p[i] = FastMath.sqrt(FastMath.pow(real, 2) + FastMath.pow(imag, 2));
	        }
		} else {
			p = new double[n];
	        for (int i = 0; i < n; i++) {
	            double real = 0;
	            double imag = 0;
	            for (int j = 0; j < n; j++) {
	                double angle = (2 * Math.PI * i * j) / n;
	                real += x[j] * FastMath.cos(angle);
	                imag += x[j] * FastMath.sin(angle);
	            }
	            p[i] = FastMath.sqrt(FastMath.pow(real, 2) + FastMath.pow(imag, 2));
	        }
		}
        return p;
	}

}
