package net.sytes.botg.array.spectrum;

import net.sytes.botg.array.math.Scalar;
import net.sytes.botg.array.math.Vec;

public class Spectrum {
	
	// Suppress default constructor for noninstantiability
	private Spectrum() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * computes the spectrogram data for {@code x} with {@code windowSize}, {@code overlapFactor} and {@code windowFunction}
	 * <br>it returns the result in form of a 2D vector, where the first dimension is time and the second the frequency
	 * @param x
	 * @param windowSize
	 * @param overlapFactor
	 * @param windowFunction
	 * @return
	 */
	public static double[][] spectrogram(double[] x, int windowSize, int overlapFactor, double[] windowFunction) {
		double[][] windows = Vec.overlapWindows(x, windowSize, overlapFactor);
		
		int w = windows.length;
				
		// apply window function if specified
		if (windowFunction != null) {
			// check if window size is correct
			if (windowFunction.length != windowSize) {
				throw new IllegalArgumentException("Window Function window must be equal to fftSampleSize");
			}
			for (int i = 0; i < w; i++) {
				windows[w] = Vec.product(windows[w], windowFunction);
			}
		}
		
		double[][] absoluteSpectrogram = new double[w][];
		
		// compute FFT for each window
		// remember min/max of spectrum
		double maxAmp = Double.MIN_VALUE;
		double minAmp = Double.MAX_VALUE;
		
		int spectrumLen = 0;
		
		for (int i = 0; i < w; i++) {
			absoluteSpectrogram[i] = singleSidedSpectrum(windows[i]);
			spectrumLen = absoluteSpectrogram[i].length;
			for (int j = 0; j < absoluteSpectrogram[i].length; j++) {
				if (absoluteSpectrogram[i][j] > maxAmp) {
					maxAmp = absoluteSpectrogram[i][j];
				}
				if (absoluteSpectrogram[i][j] < minAmp) {
					minAmp = absoluteSpectrogram[i][j];
				}
			}
		}

		boolean normalize = true;
		
		// normalization
		// avoiding divided by zero
		if (normalize) {
			double minValidAmp = 0.00000000001F;
			if (minAmp == 0) {
				minAmp = minValidAmp;
			}
	
			double[][] normalizedSpectrogram = new double[w][spectrumLen];
			
			double diff = Math.log10(maxAmp / minAmp); // perceptual difference
			for (int i = 0; i < w; i++) {
				for (int j = 0; j < spectrumLen; j++) {
					if (absoluteSpectrogram[i][j] < minValidAmp) {
						normalizedSpectrogram[i][j] = 0;
					} else {
						normalizedSpectrogram[i][j] = (Math.log10(absoluteSpectrogram[i][j] / minAmp)) / diff;
						// System.out.println(spectrogram[i][j]);
					}
				}
			}
			return normalizedSpectrogram;
		} else {
			return absoluteSpectrogram;
		}
	}
	
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
     * compute the DFT of x[] via brute force (n^2 time)
     * @param x
     * @return
     */
    public static Complex[] dft(Complex[] x) {
        int n = x.length;
        Complex ZERO = new Complex(0, 0);
        Complex[] y = new Complex[n];
        for (int k = 0; k < n; k++) {
            y[k] = ZERO;
            for (int j = 0; j < n; j++) {
                int power = (k * j) % n;
                double kth = -2 * power *  Math.PI / n;
                Complex wkj = new Complex(Math.cos(kth), Math.sin(kth));
                y[k] = y[k].add(x[j].multiply(wkj));
            }
        }
        return y;
    }
		
	/**
	 * returns the single sided spectrum p of x
	 * NOTE: input array will be trimmed to 2^n length
	 * @param x
	 * @return
	 */
	public static double[] singleSidedSpectrum(double[] x) {
		Complex[] X = fft(x, true);
		int L = X.length;	
		double[] p = new double[L / 2 + 1];
		for (int i = 0; i < L / 2 + 1; i++) {
			p[i] = X[i].abs() / L;
			if (i != 0 && i == L / 2) {
				p[i] = p[i] * 2;
			}
		}
		return p;
	}
	
	/**
	 * returns the single sided spectrum p and corresponding frequencies f of x
	 * NOTE: input array will be trimmed to 2^n length
	 * @param x
	 * @param sampleRate
	 * @return double[][] S, where S[0]: contains the frequencies and S[1] the spectrum amplitudes
	 */
	public static double[][] singleSidedSpectrum(double[] x, double sampleRate) {
		Complex[] X = fft(x, true);
		int L = X.length;		
		double[][] p = new double[2][L / 2 + 1];
		double f = 0.0;
		for (int i = 0; i < L / 2 + 1; i++) {
			p[1][i] = X[i].abs() / L;
			if (i != 0 && i < L / 2) {
				p[1][i] = p[1][i] * 2;
			}
			f = (double) i / (double) L * sampleRate;
			p[0][i] = f;			
		}
		return p;
	}
	
	public static double[][] singleSidedSpectrum2(double[] x, double sampleRate) {
		Complex[] X = fft(x, true);
		int L = X.length;		
		double[][] p = new double[2][L / 2 + 1];
		double f = 0.0;
		for (int i = 0; i < L / 2 + 1; i++) {
			p[1][i] = X[i].multiply(1.0 / (double) L).abs();
			if (i != 0 && i < L / 2) {
				p[1][i] = p[1][i] * 2;
			}
			f = (double) i / (double) L * sampleRate;
			p[0][i] = f;			
		}
		return p;
	}
	
	/**
	 * compute fft based on double input
	 * @param d
	 * @return
	 */
    public static Complex[] fft(double[] d) {
    	return fft(d, false);
    }
    
    /**
     * compute fft based on double input, where input length is trimmed to 2 * n
     * @param d
     * @param trim
     * @return
     */
    public static Complex[] fft(double[] d, boolean trim) {
    	Complex[] c;
    	if (!trim) {
    		if (!Scalar.isDivByPow(d.length, 2)) {
	            throw new IllegalArgumentException("n is not a power of 2");
	        }
    		c = toComplex(d);
    	} else {
    		int n = Scalar.nextLowerExponentForBase2(d.length);    		
    		double[] dn = new double[(int) Math.pow(2.0, (double) n)];
    		System.arraycopy(d, 0, dn, 0, dn.length);
    		c = toComplex(dn);
    	}
		return fft(c);
    }
    
    /**
     * compute the fft for specified vector in {@code d} by expanding the vector to a length of power of 2
     * @param d
     * @return
     */
    public static Complex[] fft2(double[] d) {    	
    	if (d.length < 2) {
    		throw new IllegalArgumentException("Array must be at least of length 2");
    	}
    	
    	int n = Scalar.nextGreaterExponentForBase2(d.length);
    	double[] nd = new double[(int) Math.pow(2, n)];
    	System.arraycopy(d, 0, nd, 0, d.length);
    	
    	int rem = nd.length - d.length;
    	
    	for (int i = 0; i < rem; i++) {
    		nd[i + d.length] = d[d.length - 2 - i];
    	}
    	
    	return fft(nd);
    }
	
    /**
     * Radix-2 Cooley-Tukey FFT ALGORITHM
     * compute the FFT of x[], assuming its length n is a power of 2
     * @param x
     * @see https://www.ams.org/journals/mcom/1965-19-090/S0025-5718-1965-0178586-1/
     * @return
     */
    public static Complex[] fft(Complex[] x) {
        int n = x.length;
        // base case
        if (n == 1) return new Complex[] { x[0] };

        //
        
        //if (n % 2 != 0) {
        //    throw new IllegalArgumentException("n is not a power of 2");
        //}

        // compute FFT of even terms
        Complex[] even = new Complex[n/2];
        for (int k = 0; k < n/2; k++) {
            even[k] = x[2*k];
        }
        Complex[] evenFFT = fft(even);

        // compute FFT of odd terms
        Complex[] odd  = even;  // reuse the array (to avoid n log n space)
        for (int k = 0; k < n/2; k++) {
            odd[k] = x[2*k + 1];
        }
        Complex[] oddFFT = fft(odd);

        // combine
        Complex[] y = new Complex[n];
        for (int k = 0; k < n/2; k++) {
            double kth = -2 * k * Math.PI / n;
            Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
            y[k]       = evenFFT[k].add(wk.multiply(oddFFT[k]));
            y[k + n/2] = evenFFT[k].subtract(wk.multiply(oddFFT[k]));
        }
        return y;
    }
    
	/**
	 * convert double array to complex array
	 * @param d
	 * @return
	 */
    public static Complex[] toComplex(double[] d) {
		Complex[] c = new Complex[d.length];
		//for (int i = 0; i < c.length / 2; i++) {
		for (int i = 0; i < c.length; i++) {
			//c[2 * i + 1] = new Complex(d[2 * i + 1], 0.0);
			//c[2 * i] = new Complex(d[2 * i], 0.0);
			c[i] = new Complex(d[i], 0.0);
		}
    	return c;
	}

	/**
	 * compute the inverse FFT of x[], assuming its length n is a power of 2
	 * @param x
	 * @return
	 */
    public static Complex[] ifft(Complex[] x) {
        int n = x.length;
        Complex[] y = new Complex[n];

        // take conjugate
        for (int i = 0; i < n; i++) {
            y[i] = x[i].conjugate();
        }

        // compute forward FFT
        y = fft(y);

        // take conjugate again
        for (int i = 0; i < n; i++) {
            y[i] = y[i].conjugate();
        }

        // divide by n
        for (int i = 0; i < n; i++) {
            y[i] = y[i].multiply(1.0 / n);
        }

        return y;

    }

    /**
     * compute the circular convolution of x and y
     * @param x
     * @param y
     * @return
     */
    public static Complex[] cconvolve(Complex[] x, Complex[] y) {

        // should probably pad x and y with 0s so that they have same length
        // and are powers of 2
        if (x.length != y.length) {
            throw new IllegalArgumentException("Dimensions don't agree");
        }

        int n = x.length;

        // compute FFT of each sequence
        Complex[] a = fft(x);
        Complex[] b = fft(y);

        // point-wise multiply
        Complex[] c = new Complex[n];
        for (int i = 0; i < n; i++) {
            c[i] = a[i].multiply(b[i]);
        }

        // compute inverse FFT
        return ifft(c);
    }

    /**
     * compute the linear convolution of x and y
     * @param x
     * @param y
     * @return
     */
    public static Complex[] convolve(Complex[] x, Complex[] y) {
        Complex ZERO = new Complex(0, 0);

        Complex[] a = new Complex[2*x.length];
        for (int i = 0;        i <   x.length; i++) a[i] = x[i];
        for (int i = x.length; i < 2*x.length; i++) a[i] = ZERO;

        Complex[] b = new Complex[2*y.length];
        for (int i = 0;        i <   y.length; i++) b[i] = y[i];
        for (int i = y.length; i < 2*y.length; i++) b[i] = ZERO;

        return cconvolve(a, b);
    }
}
