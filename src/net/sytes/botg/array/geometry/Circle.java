package net.sytes.botg.array.geometry;

import net.sytes.botg.array.math.Vec;

public class Circle extends Curve2D {

	private double D;
	private double x_0;
	private double y_0;

	public Circle() {
		this(new Builder());
	}

	private Circle(Builder builder) {
		
		this.D = builder.D;
		this.x_0 = builder.x_0;
		this.y_0 = builder.y_0;
		
		// create function handles
		this.xHandle = new FunctionHandle() {
			@Override
			public double apply(double t1) {
				double theta = t1 * 2 * Math.PI;
				return x_0 + D / 2.0 * Math.cos(theta);
			}
		};
		
		this.yHandle = new FunctionHandle() {
			@Override
			public double apply(double t1) {
				double theta = t1 * 2 * Math.PI;
				return y_0 + D / 2.0 * Math.sin(theta);
			}
		};
		
	}

	public static class Builder {
		private double D;
		private double x_0;
		private double y_0;

		public Builder D(double D) {
			this.D = D;
			return this;
		}

		public Builder x_0(double x_0) {
			this.x_0 = x_0;
			return this;
		}

		public Builder y_0(double y_0) {
			this.y_0 = y_0;
			return this;
		}

		public Circle build() {
			return new Circle(this);
		}
	}
	
	@Override
	public void create(double[] t) {
		// in order to close the circle line, the entries at first index, are also the last entry, which is added additionally
		this.t = new double[t.length + 1];
		System.arraycopy(t, 0, this.t, 0, t.length);
		this.t[t.length] = t[0];
		this.x = new double[t.length + 1];
		this.y = new double[t.length + 1];
		int i = 0;
		for (i = 0; i < t.length; i++) {
			this.x[i] = this.xHandle.apply(this.t[i]);
			this.y[i] = this.yHandle.apply(this.t[i]);
		}
		this.x[i] = this.x[0];
		this.y[i] = this.y[0];
	}
	
	public static double[][] circle(double D, double x_0, double y_0, int res){
		double[][] xyData = new double[2][];
		
		Circle circle = new Circle.Builder()
				.D(D)
				.x_0(x_0)
				.y_0(y_0)
				.build();
		
		double[] t = Vec.linspace(0.0, 1.0, res);
		
		circle.create(t);
		
		xyData[0] = circle.x();
		xyData[1] = circle.y();
		
		return xyData;
	}
	
}
