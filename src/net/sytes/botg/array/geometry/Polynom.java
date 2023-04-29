package net.sytes.botg.array.geometry;

public class Polynom extends Curve2D {

	private double[] weights;
	
	public Polynom() {
		this(new Builder());
	}	
	
	private Polynom(Builder builder) {
		this.weights = builder.weights;
		
		// create function handles
		this.xHandle = new FunctionHandle() {
			@Override
			public double apply(double t1) {
				return t1;
			}
		};
		
		this.yHandle = new FunctionHandle() {
			@Override
			public double apply(double t1) {
				double p = 0;
				for (int i = 0; i < weights.length; i++) {
					p = p + weights[i] * Math.pow(t1, i);
				}
				return p;
			}
		};
		
	}
	
	public static class Builder {
		
		private double[] weights = new double[2];
				
		public Builder weights(double[] weights) {
			this.weights = weights;
			return this;
		}

		public Polynom build() {
			return new Polynom(this);
		}
		
	}
	
	public int degree() {
		return this.weights.length - 1;
	}	
	
}
