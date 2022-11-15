package net.sytes.botg.array.spectrum;

import org.apache.commons.math3.complex.Complex;

import net.sytes.botg.array.math.Scalar;

/******************************************************************************
 * Copyright © 2000–2019, Robert Sedgewick and Kevin Wayne.
 *	Last updated: Tue Jan 14 09:42:25 EST 2020. 
 * 
 * Modified 2021, Jonas Hillenbrand
 * 
 *  Compilation:  javac FFT.java
 *  Execution:    java FFT n
 *  Dependencies: Complex.java
 *
 *  Compute the FFT and inverse FFT of a length n complex sequence
 *  using the radix 2 Cooley-Tukey algorithm.

 *  Bare bones implementation that runs in O(n log n) time and O(n)
 *  space. Our goal is to optimize the clarity of the code, rather
 *  than performance.
 *
 *  This implementation uses the primitive root of unity w = e^(-2 pi i / n).
 *  Some resources use w = e^(2 pi i / n).
 *
 *  Reference: https://www.cs.princeton.edu/~wayne/kleinberg-tardos/pdf/05DivideAndConquerII.pdf
 *
 *  Limitations
 *  -----------
 *   -  assumes n is a power of 2
 *
 *   -  not the most memory efficient algorithm (because it uses
 *      an object type for representing complex numbers and because
 *      it re-allocates memory for the subarray, instead of doing
 *      in-place or reusing a single temporary array)
 *  
 *  For an in-place radix 2 Cooley-Tukey FFT, see
 *  https://introcs.cs.princeton.edu/java/97data/InplaceFFT.java.html
 *
 ******************************************************************************/

public class FFT {

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
	 * @return
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
    		int n = Scalar.closestExponentForBase2(d.length);    		
    		double[] dn = new double[(int) Math.pow(2.0, (double) n)];
    		System.arraycopy(d, 0, dn, 0, dn.length);
    		c = toComplex(dn);
    	}
		return fft(c);
    }
	
    /**
     *  Radix-2 Cooley-Tukey FFT ALGOTRIHM
     *  compute the FFT of x[], assuming its length n is a power of 2
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

	// compute the inverse FFT of x[], assuming its length n is a power of 2
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

    // compute the circular convolution of x and y
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


    // compute the linear convolution of x and y
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

    // compute the DFT of x[] via brute force (n^2 time)
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
}