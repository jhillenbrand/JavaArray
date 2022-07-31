package net.sytes.botg.array.spectrum;

public class Complex {

    private final double imaginary;
    private final double real;
    
    /** Mask used to clear the non-sign part of a long. */
    private static final long MASK_NON_SIGN_LONG = 0x7fffffffffffffffl;
    
    /**
     * instantiate a complex number by specifying only the real part.
     *
     * @param real real part.
     */
    public Complex(double real) {
        this(real, 0.0);
    }

    /**
     * instantiate a complex number by specifying the real and imaginary parts.
     *
     * @param real real part.
     * @param imaginary imaginary part.
     */
    public Complex(double real, double imaginary) {
        this.real = real;
        this.imaginary = imaginary;
    }
    
    public double abs() {
    	Double.longBitsToDouble(MASK_NON_SIGN_LONG & Double.doubleToRawLongBits(x));
    }
	
}
