package net.sytes.botg.array.spectrum;

public class Complex {

    private final double imaginary;
    private final double real;
    
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
    	return Math.sqrt(this.real * this.real + this.imaginary * this.imaginary);
    }
    
    public Complex multiply(Complex c) {
    	return new Complex(this.real * c.real - this.imaginary * c.imaginary, this.real * c.imaginary + this.imaginary * c.real);
    }
    
    public Complex multiply(double d) {
    	return this.multiply(new Complex(d));
    }
    
    public Complex add(Complex c) {
    	return new Complex(this.real + c.real, this.imaginary + c.imaginary);
    }
    
    public Complex subtract(Complex c) {
    	return new Complex(this.real - c.real, this.imaginary - c.imaginary);
    }
    
    public Complex square() {
    	return new Complex(Math.pow(this.real, 2) - Math.pow(this.imaginary, 2), 2 * this.real * this.imaginary);
    }
	
    public Complex conjugate() {
    	return new Complex(this.real, -this.imaginary);
    }
    
    public double getImaginary() {
		return this.imaginary;
	}

	public double getReal() {
		return this.real;
	}

	@Override
    public String toString() {
    	if (this.imaginary < 0) {
        	return "" + this.real + " - " + Math.abs(this.imaginary) + "i";
    	} else {
        	return "" + this.real + " + " + this.imaginary + "i";
    	}
    }
    
}
