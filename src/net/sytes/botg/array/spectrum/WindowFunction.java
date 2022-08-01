package net.sytes.botg.array.spectrum;

import net.sytes.botg.array.ArUtils;

/**
 * Utility class that provides various window types that can be applied on sample vectors
 * See <a href="https://en.wikipedia.org/wiki/Window_function">Wikipedia</a> for more info
 * @author jhillenb
 *
 */
public class WindowFunction {

	// Suppress default constructor for noninstantiability
	private WindowFunction() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	public enum WindowType {
		RECTANGULAR, TRIANGULAR, PARZEN, WELCH, SINE, POWER_OF_SINE, HANN, HAMMING, BLACKMAN,
		NUTALL, BLACKMAN_NUTALL, BLACKMAN_HARRIS, FLAT_TOP, RIFE_VINCENT, GAUSSIAN, GENERALIZED_NORMAL,
		TUKEY, PLANCK_TAPER, SLEPIAN, KAISER, DOLPH_CHEBYSHEV, ULTRASPHERICAL, POISSON, BARTLET_HANN,
		PLANCK_BESSEL, GAP, LANCZOS
	}
	
	/**
	 * generates a single window weight of the {@code WindowType} {@code wt} for the specified index {@code i} in window of length {@code n}
	 * <br>
	 * <br>returns -1.0 if {@code WindowType} was not defined yet
	 * @param wt
	 * @param i
	 * @param n
	 * @return
	 */
	public double generate(WindowType wt, int i, int n) {
		switch (wt) {
			case BARTLET_HANN:
				break;
			case BLACKMAN:
				break;
			case BLACKMAN_HARRIS:
				return blackmanHarris(i, n);
				
			case BLACKMAN_NUTALL:
				break;
			case DOLPH_CHEBYSHEV:
				break;
			case FLAT_TOP:
				break;
			case GAP:
				break;
			case GAUSSIAN:
				break;
			case GENERALIZED_NORMAL:
				break;
			case HAMMING:
				return hamming(i, n);
				
			case HANN:
				return hann(i, n);
				
			case KAISER:
				break;
			case LANCZOS:
				break;
			case NUTALL:
				break;
			case PARZEN:
				break;
			case PLANCK_BESSEL:
				break;
			case PLANCK_TAPER:
				break;
			case POISSON:
				break;
			case POWER_OF_SINE:
				break;
			case RECTANGULAR:
				return rectangular(i, n);
				
			case RIFE_VINCENT:
				break;
			case SINE:
				break;
			case SLEPIAN:
				break;
			case TRIANGULAR:
				return triangular(i, n);
				
			case TUKEY:
				break;
			case ULTRASPHERICAL:
				break;
			case WELCH:
				return welch(i, n);
		}
		return -1.0;
	}
	
	/**
	 * generates a window of size {@code n} for {@code WindowType} {@code wt}
	 * <br>
	 * <br>returns {@code null} if {@code wt} was not defined yet
	 * @param wt
	 * @param n
	 * @return
	 */
	public double[] generate(WindowType wt, int n) {
		switch (wt) {
			case BARTLET_HANN:
				break;
			case BLACKMAN:
				break;
			case BLACKMAN_HARRIS:
				return blackmanHarris(n);
				
			case BLACKMAN_NUTALL:
				break;
			case DOLPH_CHEBYSHEV:
				break;
			case FLAT_TOP:
				break;
			case GAP:
				break;
			case GAUSSIAN:
				break;
			case GENERALIZED_NORMAL:
				break;
			case HAMMING:
				return hamming(n);
				
			case HANN:
				return hann(n);
				
			case KAISER:
				break;
			case LANCZOS:
				break;
			case NUTALL:
				break;
			case PARZEN:
				break;
			case PLANCK_BESSEL:
				break;
			case PLANCK_TAPER:
				break;
			case POISSON:
				break;
			case POWER_OF_SINE:
				break;
			case RECTANGULAR:
				return rectangular(n);
				
			case RIFE_VINCENT:
				break;
			case SINE:
				break;
			case SLEPIAN:
				break;
			case TRIANGULAR:
				return triangular(n);
				
			case TUKEY:
				break;
			case ULTRASPHERICAL:
				break;
			case WELCH:
				return welch(n);		
		}
		return null;
	}
	
	/**
	 * applies the {@code WindowType} {@code wt} on the input vector {@code x} and<br>returns the weighted result as {@code double[]}
	 * <br>
	 * <br>returns {@code null} if {@code wt} was not defined yet 
	 * @param wt
	 * @param x
	 * @return
	 */
	public double[] apply(WindowType wt, double[] x) {
		int n = x.length;
		double[] x_w = new double[n];
		switch (wt) {
			case BARTLET_HANN:
				break;
			case BLACKMAN:
				break;
			case BLACKMAN_HARRIS:
				for (int i = 0; i < n; i++) {
					x_w[i] = x[i] * blackmanHarris(i, n);
				}
				return x_w;
				
			case BLACKMAN_NUTALL:
				break;
			case DOLPH_CHEBYSHEV:
				break;
			case FLAT_TOP:
				break;
			case GAP:
				break;
			case GAUSSIAN:
				break;
			case GENERALIZED_NORMAL:
				break;
			case HAMMING:
				for (int i = 0; i < n; i++) {
					x_w[i] = x[i] * hamming(i, n);
				}
				return x_w;
				
			case HANN:
				for (int i = 0; i < n; i++) {
					x_w[i] = x[i] * hann(i, n);
				}
				return x_w;
				
			case KAISER:
				break;
			case LANCZOS:
				break;
			case NUTALL:
				break;
			case PARZEN:
				break;
			case PLANCK_BESSEL:
				break;
			case PLANCK_TAPER:
				break;
			case POISSON:
				break;
			case POWER_OF_SINE:
				break;
			case RECTANGULAR:
				for (int i = 0; i < n; i++) {
					x_w[i] = x[i] * rectangular(i, n);
				}
				return x_w;
				
			case RIFE_VINCENT:
				break;
			case SINE:
				break;
			case SLEPIAN:
				break;
			case TRIANGULAR:
				for (int i = 0; i < n; i++) {
					x_w[i] = x[i] * triangular(i, n);
				}
				return x_w;
				
			case TUKEY:
				break;
			case ULTRASPHERICAL:
				break;
			case WELCH:
				for (int i = 0; i < n; i++) {
					x_w[i] = x[i] * welch(i, n);
				}
				return x_w;	
		}
		return null;
	}
	
	/**
	 * applies the {@code WindowType} {@code wt} for a single sample {@code x_i} at index {@code i} of a window of length {@code n}
	 * @param wt
	 * @param i
	 * @param n
	 * @param x_i
	 * @return
	 */
	public double apply(WindowType wt, int i, int n, double x_i) {
		switch (wt) {
			case BARTLET_HANN:
				break;
			case BLACKMAN:
				break;
			case BLACKMAN_HARRIS:
				return x_i * blackmanHarris(i, n);
				
			case BLACKMAN_NUTALL:
				break;
			case DOLPH_CHEBYSHEV:
				break;
			case FLAT_TOP:
				break;
			case GAP:
				break;
			case GAUSSIAN:
				break;
			case GENERALIZED_NORMAL:
				break;
			case HAMMING:
				return x_i * hamming(i, n);
				
			case HANN:
				return x_i * hann(i, n);
				
			case KAISER:
				break;
			case LANCZOS:
				break;
			case NUTALL:
				break;
			case PARZEN:
				break;
			case PLANCK_BESSEL:
				break;
			case PLANCK_TAPER:
				break;
			case POISSON:
				break;
			case POWER_OF_SINE:
				break;
			case RECTANGULAR:
				return x_i * rectangular(i, n);
				
			case RIFE_VINCENT:
				break;
			case SINE:
				break;
			case SLEPIAN:
				break;
			case TRIANGULAR:
				return x_i * triangular(i, n);
				
			case TUKEY:
				break;
			case ULTRASPHERICAL:
				break;
			case WELCH:
				return x_i * welch(i, n);
		}
		return -1.0;
	}
	
	private static double rectangular(int i, int n) {
		return 1.0;
	}
	
	private static double[] rectangular(int n) {
		return ArUtils.ones(n);
	}
	
	private static double triangular(int i, int n) {
		return 1.0 - Math.abs((i - n / 2.0) / (n / 2.0));
	}
	
	private static double[] triangular(int n) {
		double[] w = new double[n];
		for (int i = 0; i < n; i++) {
			w[i] = triangular(i, n);
		}
		return w;
	}
	
	private static double hann(int i, int n) {
		return 0.5 * (1 - Math.cos(2 * Math.PI * i / n)); 
	}
	
	private static double[] hann(int n) {
		double[] w = new double[n];
		for (int i = 0; i < n; i++) {
			w[i] = hann(i, n);
		}
		return w;
	}
	
	private static double hamming(int i, int n) {
		return 0.53836 - (1 - 0.46164) * Math.cos(2 * Math.PI * i / n);
	}
	
	private static double[] hamming(int n) {
		double[] w = new double[n];
		for (int i = 0; i < n; i++) {
			w[i] = hamming(i, n);
		}
		return w;
	}
	
	private static double blackmanHarris(int i, int n) {
		return 0.35875 - 0.48829 * Math.cos(2 * Math.PI * i / n) + 0.14128 * Math.cos(4 * Math.PI * i / n) - 0.01168 * Math.cos(6 * Math.PI * i / n); 
	}
	
	private static double[] blackmanHarris(int n) {
		double[] w = new double[n];
		for (int i = 0; i < n; i++) {
			w[i] = blackmanHarris(i, n);
		}
		return w;
	}
	
	private static double welch(int i, int n) {
		return 1 - Math.pow((i - n / 2.0) / (n / 2.0), 2);
	}
	
	private static double[] welch(int n) {
		double[] w = new double[n];
		for (int i = 0; i < n; i++) {
			w[i] = welch(i, n);
		}
		return w;
	}
	
}
