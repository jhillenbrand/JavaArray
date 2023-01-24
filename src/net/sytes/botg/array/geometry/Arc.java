package net.sytes.botg.array.geometry;

public class Arc extends Curve2D {

	private double xs = 0.0;
	private double ys = 1.0;
	private double r = 1.0;	
	private double xe = 1.0;
	private double ye = 0.0;
	
	private double xc = 0.0;
	private double yc = 0.0;
	
	private double theta1 = 0.0;
	private double theta2 = 0.0;

	public Arc() {
		this(new Builder());
	}
		
	private Arc(Builder builder) {
		this.r = builder.r;
		this.xc = builder.xc;
		this.yc = builder.yc;
		this.xs = builder.x1;
		this.ys = builder.y1;
		this.xe = builder.x2;
		this.ye = builder.y2;
		
		// compute start and stop angles theta0 and theta1
		// differentiate quadrants
		// Q1 -> 0°-90°
		if(this.xs >= 0 && this.ys >= 0) {
			this.theta1 = Math.atan((this.xs - this.xc) / (this.ys - this.yc)); 
		}
		if (this.xe >= 0 && this.ye >= 0) {
			this.theta2 = Math.atan((this.xe - this.xc) / (this.ye - this.yc));
		}
		// Q2 -> 90°-180°
		if (this.xs >= 0 && this.ys < 0) {
			this.theta1 = Math.PI / 2 + Math.atan((this.ys - this.yc) / (this.xs - this.xc));
		}
		if (this.xe >= 0 && this.ye < 0) {
			this.theta2 = Math.PI / 2 + Math.atan((this.ye - this.yc) / (this.xe - this.xc));
		}
		// Q3 -> 180°-270°
		if (this.xs < 0 && this.ys < 0) {
			this.theta1 = Math.PI + Math.atan((this.xc - this.xs) / (this.yc - this.ys));
		}
		if (this.xe < 0 && this.ye < 0) {
			this.theta2 = Math.PI + Math.atan((this.xc - this.xe) / (this.yc - this.ye));
		}
		// Q4 -> 270°-360° 
		if (this.xs < 0 && this.ys >= 0) {
			this.theta1 = 3 / 2 * Math.PI / 2 + Math.atan((this.ys - this.yc) / (this.xc - this.xs));
		}
		if (this.xe < 0 && this.ye >= 0) {
			this.theta2 = 3 / 2 * Math.PI / 2 + Math.atan((this.ye - this.yc) / (this.xc - this.xe));
		}
		
		// create function handles
		this.xHandle = new FunctionHandle() {
			@Override
			public double apply(double t1) {
				double theta = theta1 + t1 * (theta2 - theta1);
				return xc + r * Math.cos(theta);
			}
		};
		
		this.yHandle = new FunctionHandle() {
			@Override
			public double apply(double t1) {
				double theta = t1 * 2 * Math.PI;
				return yc + r * Math.sin(theta);
			}
		};
		
	}
	
	public static class Builder {
		private double r;
		private double xc;
		private double yc;
		private double x1;
		private double y1;
		private double x2;
		private double y2;

		public Builder r(double r) {
			this.r = r;
			return this;
		}

		public Builder xc(double xc) {
			this.xc = xc;
			return this;
		}

		public Builder yc(double yc) {
			this.yc = yc;
			return this;
		}

		public Builder x1(double x1) {
			this.x1 = x1;
			return this;
		}

		public Builder y1(double y1) {
			this.y1 = y1;
			return this;
		}

		public Builder x2(double x2) {
			this.x2 = x2;
			return this;
		}

		public Builder y2(double y2) {
			this.y2 = y2;
			return this;
		}

		public Arc build() {
			return new Arc(this);
		}
	}
	
}
