package net.sytes.botg.array.spectrum;

import net.sytes.botg.array.math.Mat;
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
    
    /**
     * Algorithm for Singular Spectrum Analysis based on <a href="https://valeman.medium.com/singular-spectrum-analysis-a-hidden-treasure-of-time-series-forecasting-b8ae8b30a948">Link</a>
     * @param y
     * @param w
     * @return
     */
    public static double[][] singularSpectrumAnalysis(double[] x, int L){
    	
    	// 1) embedding step    	
    	double[][] X = Mat.trajectoryMatrix(x, L);
    	
    	// 2) singular value decomposition (SVD)
        // Derived from LINPACK code / Jama
        // Initialize.
        int m = X.length;
        int n = X[0].length;

        /* Apparently the failing cases are only a proper subset of (m<n), 
	  	 so let's not throw error.  Correct fix to come later?
	        if (m<n) {
	  	  throw new IllegalArgumentException("Jama SVD only works for m >= n"); }
	        */
        int nu = Math.min(m,n);
        double[] s = new double [Math.min(m+1,n)];
        double[][] U = new double [m][nu];
        double[][] V = new double [n][n];
        double[] e = new double [n];
        double[] work = new double [m];
        boolean wantu = true;
        boolean wantv = true;

        // Reduce A to bidiagonal form, storing the diagonal elements
        // in s and the super-diagonal elements in e.

        int nct = Math.min(m-1,n);
        int nrt = Math.max(0,Math.min(n-2,m));
        for (int k = 0; k < Math.max(nct,nrt); k++) {
           if (k < nct) {

              // Compute the transformation for the k-th column and
              // place the k-th diagonal in s[k].
              // Compute 2-norm of k-th column without under/overflow.
              s[k] = 0;
              for (int i = k; i < m; i++) {
                 s[k] = Scalar.hypot(s[k],X[i][k]);
              }
              if (s[k] != 0.0) {
                 if (X[k][k] < 0.0) {
                    s[k] = -s[k];
                 }
                 for (int i = k; i < m; i++) {
                    X[i][k] /= s[k];
                 }
                 X[k][k] += 1.0;
              }
              s[k] = -s[k];
           }
           for (int j = k+1; j < n; j++) {
              if ((k < nct) & (s[k] != 0.0))  {

              // Apply the transformation.

                 double t = 0;
                 for (int i = k; i < m; i++) {
                    t += X[i][k]*X[i][j];
                 }
                 t = -t/X[k][k];
                 for (int i = k; i < m; i++) {
                    X[i][j] += t*X[i][k];
                 }
              }

              // Place the k-th row of X into e for the
              // subsequent calculation of the row transformation.

              e[j] = X[k][j];
           }
           if (wantu & (k < nct)) {

              // Place the transformation in U for subsequent back
              // multiplication.

              for (int i = k; i < m; i++) {
                 U[i][k] = X[i][k];
              }
           }
           if (k < nrt) {

              // Compute the k-th row transformation and place the
              // k-th super-diagonal in e[k].
              // Compute 2-norm without under/overflow.
              e[k] = 0;
              for (int i = k+1; i < n; i++) {
                 e[k] = Scalar.hypot(e[k],e[i]);
              }
              if (e[k] != 0.0) {
                 if (e[k+1] < 0.0) {
                    e[k] = -e[k];
                 }
                 for (int i = k+1; i < n; i++) {
                    e[i] /= e[k];
                 }
                 e[k+1] += 1.0;
              }
              e[k] = -e[k];
              if ((k+1 < m) & (e[k] != 0.0)) {

              // Apply the transformation.

                 for (int i = k+1; i < m; i++) {
                    work[i] = 0.0;
                 }
                 for (int j = k+1; j < n; j++) {
                    for (int i = k+1; i < m; i++) {
                       work[i] += e[j]*X[i][j];
                    }
                 }
                 for (int j = k+1; j < n; j++) {
                    double t = -e[j]/e[k+1];
                    for (int i = k+1; i < m; i++) {
                       X[i][j] += t*work[i];
                    }
                 }
              }
              if (wantv) {

              // Place the transformation in V for subsequent
              // back multiplication.

                 for (int i = k+1; i < n; i++) {
                    V[i][k] = e[i];
                 }
              }
           }
        }

        // Set up the final bidiagonal matrix or order p.

        int p = Math.min(n,m+1);
        if (nct < n) {
           s[nct] = X[nct][nct];
        }
        if (m < p) {
           s[p-1] = 0.0;
        }
        if (nrt+1 < p) {
           e[nrt] = X[nrt][p-1];
        }
        e[p-1] = 0.0;

        // If required, generate U.

        if (wantu) {
           for (int j = nct; j < nu; j++) {
              for (int i = 0; i < m; i++) {
                 U[i][j] = 0.0;
              }
              U[j][j] = 1.0;
           }
           for (int k = nct-1; k >= 0; k--) {
              if (s[k] != 0.0) {
                 for (int j = k+1; j < nu; j++) {
                    double t = 0;
                    for (int i = k; i < m; i++) {
                       t += U[i][k]*U[i][j];
                    }
                    t = -t/U[k][k];
                    for (int i = k; i < m; i++) {
                       U[i][j] += t*U[i][k];
                    }
                 }
                 for (int i = k; i < m; i++ ) {
                    U[i][k] = -U[i][k];
                 }
                 U[k][k] = 1.0 + U[k][k];
                 for (int i = 0; i < k-1; i++) {
                    U[i][k] = 0.0;
                 }
              } else {
                 for (int i = 0; i < m; i++) {
                    U[i][k] = 0.0;
                 }
                 U[k][k] = 1.0;
              }
           }
        }

        // If required, generate V.

        if (wantv) {
           for (int k = n-1; k >= 0; k--) {
              if ((k < nrt) & (e[k] != 0.0)) {
                 for (int j = k+1; j < nu; j++) {
                    double t = 0;
                    for (int i = k+1; i < n; i++) {
                       t += V[i][k]*V[i][j];
                    }
                    t = -t/V[k+1][k];
                    for (int i = k+1; i < n; i++) {
                       V[i][j] += t*V[i][k];
                    }
                 }
              }
              for (int i = 0; i < n; i++) {
                 V[i][k] = 0.0;
              }
              V[k][k] = 1.0;
           }
        }

        // Main iteration loop for the singular values.

        int pp = p-1;
        int iter = 0;
        double eps = Math.pow(2.0,-52.0);
        double tiny = Math.pow(2.0,-966.0);
        while (p > 0) {
           int k,kase;

           // Here is where a test for too many iterations would go.

           // This section of the program inspects for
           // negligible elements in the s and e arrays.  On
           // completion the variables kase and k are set as follows.

           // kase = 1     if s(p) and e[k-1] are negligible and k<p
           // kase = 2     if s(k) is negligible and k<p
           // kase = 3     if e[k-1] is negligible, k<p, and
           //              s(k), ..., s(p) are not negligible (qr step).
           // kase = 4     if e(p-1) is negligible (convergence).

           for (k = p-2; k >= -1; k--) {
              if (k == -1) {
                 break;
              }
              if (Math.abs(e[k]) <=
                    tiny + eps*(Math.abs(s[k]) + Math.abs(s[k+1]))) {
                 e[k] = 0.0;
                 break;
              }
           }
           if (k == p-2) {
              kase = 4;
           } else {
              int ks;
              for (ks = p-1; ks >= k; ks--) {
                 if (ks == k) {
                    break;
                 }
                 double t = (ks != p ? Math.abs(e[ks]) : 0.) + 
                            (ks != k+1 ? Math.abs(e[ks-1]) : 0.);
                 if (Math.abs(s[ks]) <= tiny + eps*t)  {
                    s[ks] = 0.0;
                    break;
                 }
              }
              if (ks == k) {
                 kase = 3;
              } else if (ks == p-1) {
                 kase = 1;
              } else {
                 kase = 2;
                 k = ks;
              }
           }
           k++;

           // Perform the task indicated by kase.

           switch (kase) {

              // Deflate negligible s(p).

              case 1: {
                 double f = e[p-2];
                 e[p-2] = 0.0;
                 for (int j = p-2; j >= k; j--) {
                    double t = Scalar.hypot(s[j],f);
                    double cs = s[j]/t;
                    double sn = f/t;
                    s[j] = t;
                    if (j != k) {
                       f = -sn*e[j-1];
                       e[j-1] = cs*e[j-1];
                    }
                    if (wantv) {
                       for (int i = 0; i < n; i++) {
                          t = cs*V[i][j] + sn*V[i][p-1];
                          V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
                          V[i][j] = t;
                       }
                    }
                 }
              }
              break;

              // Split at negligible s(k).

              case 2: {
                 double f = e[k-1];
                 e[k-1] = 0.0;
                 for (int j = k; j < p; j++) {
                    double t = Scalar.hypot(s[j],f);
                    double cs = s[j]/t;
                    double sn = f/t;
                    s[j] = t;
                    f = -sn*e[j];
                    e[j] = cs*e[j];
                    if (wantu) {
                       for (int i = 0; i < m; i++) {
                          t = cs*U[i][j] + sn*U[i][k-1];
                          U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
                          U[i][j] = t;
                       }
                    }
                 }
              }
              break;

              // Perform one qr step.

              case 3: {

                 // Calculate the shift.
     
                 double scale = Math.max(Math.max(Math.max(Math.max(
                         Math.abs(s[p-1]),Math.abs(s[p-2])),Math.abs(e[p-2])), 
                         Math.abs(s[k])),Math.abs(e[k]));
                 double sp = s[p-1]/scale;
                 double spm1 = s[p-2]/scale;
                 double epm1 = e[p-2]/scale;
                 double sk = s[k]/scale;
                 double ek = e[k]/scale;
                 double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
                 double c = (sp*epm1)*(sp*epm1);
                 double shift = 0.0;
                 if ((b != 0.0) | (c != 0.0)) {
                    shift = Math.sqrt(b*b + c);
                    if (b < 0.0) {
                       shift = -shift;
                    }
                    shift = c/(b + shift);
                 }
                 double f = (sk + sp)*(sk - sp) + shift;
                 double g = sk*ek;
     
                 // Chase zeros.
     
                 for (int j = k; j < p-1; j++) {
                    double t = Scalar.hypot(f,g);
                    double cs = f/t;
                    double sn = g/t;
                    if (j != k) {
                       e[j-1] = t;
                    }
                    f = cs*s[j] + sn*e[j];
                    e[j] = cs*e[j] - sn*s[j];
                    g = sn*s[j+1];
                    s[j+1] = cs*s[j+1];
                    if (wantv) {
                       for (int i = 0; i < n; i++) {
                          t = cs*V[i][j] + sn*V[i][j+1];
                          V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
                          V[i][j] = t;
                       }
                    }
                    t = Scalar.hypot(f,g);
                    cs = f/t;
                    sn = g/t;
                    s[j] = t;
                    f = cs*e[j] + sn*s[j+1];
                    s[j+1] = -sn*e[j] + cs*s[j+1];
                    g = sn*e[j+1];
                    e[j+1] = cs*e[j+1];
                    if (wantu && (j < m-1)) {
                       for (int i = 0; i < m; i++) {
                          t = cs*U[i][j] + sn*U[i][j+1];
                          U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
                          U[i][j] = t;
                       }
                    }
                 }
                 e[p-2] = f;
                 iter = iter + 1;
              }
              break;

              // Convergence.

              case 4: {

                 // Make the singular values positive.
     
                 if (s[k] <= 0.0) {
                    s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
                    if (wantv) {
                       for (int i = 0; i <= pp; i++) {
                          V[i][k] = -V[i][k];
                       }
                    }
                 }
     
                 // Order the singular values.
     
                 while (k < pp) {
                    if (s[k] >= s[k+1]) {
                       break;
                    }
                    double t = s[k];
                    s[k] = s[k+1];
                    s[k+1] = t;
                    if (wantv && (k < n-1)) {
                       for (int i = 0; i < n; i++) {
                          t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
                       }
                    }
                    if (wantu && (k < m-1)) {
                       for (int i = 0; i < m; i++) {
                          t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
                       }
                    }
                    k++;
                 }
                 iter = 0;
                 p--;
              }
              break;
           }
        }
    	
    	return null;
    }
    
}
