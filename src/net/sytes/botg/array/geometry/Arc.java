package net.sytes.botg.array.geometry;

import net.sytes.botg.array.math.Scalar;

/**
 * generates the coordinate points of an arc defined by 3 points (center, start and end)
 * <img src="Arc.png"></img>
 * @author hillenbrand
 *
 */
public class Arc extends Curve2D {

	private double xs;
	private double ys;
	private double xe;
	private double ye;	
	private double xc;
	private double yc;

	private double r = 0.0;	
	private double theta1 = 0.0;
	private double theta2 = 0.0;

	private static final double R_TOL = 1E-9;
	
	public Arc() {
		this(new Builder());
	}
		
	private Arc(Builder builder) {
		
		this.xc = builder.xc;
		this.yc = builder.yc;
		this.xs = builder.xs;
		this.ys = builder.ys;
		this.xe = builder.xe;
		this.ye = builder.ye;
		
		this.r = Scalar.norm(this.xc - this.xs, this.yc - this.ys);
		double r2 = Scalar.norm(this.xc - this.xe, this.yc - this.ye);
		if (Math.abs(r2 - this.r) > R_TOL) {
			throw new IllegalArgumentException("The starting and end point are not the same distance from the specified center (R_TOL=" + R_TOL + "). Make sure you are specifying a valid " + this.getClass().getSimpleName() + "!");
		}
		
		// compute start and stop angles theta0 and theta1
		// differentiate quadrants
		// Q1 -> 0°-90°
		if(this.xs >= this.xc && this.ys >= this.yc) {
			if (this.xs - this.xc == 0) {
				this.theta1 =  Math.PI / 2.0;
			} else {
				this.theta1 = Math.atan((this.ys - this.yc) / (this.xs - this.xc)); 
			}
		}
		if (this.xe >= this.xc && this.ye >= this.yc) {
			if (this.xe - this.xc == 0) {
				this.theta2 = Math.PI / 2.0;
			} else {
				this.theta2 = Math.atan((this.ye - this.yc) / (this.xe - this.xc));
			}
		}
		// Q2 -> 90°-180°
		if (this.xs < this.xc && this.ys >= this.yc) {
			if (this.ys - this.yc == 0){
				this.theta1 = Math.PI;
			} else {
				this.theta1 = Math.PI / 2.0 + Math.atan((this.xc - this.xs) / (this.ys - this.yc));
			}
		}
		if (this.xe < this.xc && this.ye >= this.yc) {
			if (this.ye - this.yc == 0) {
				this.theta2 = Math.PI;
			} else {
				this.theta2 = Math.PI / 2.0 + Math.atan((this.xc - this.xe) / (this.ye - this.yc));				
			}
		}
		// Q3 -> 180°-270°
		if (this.xs < this.xc && this.ys < this.yc) {
			if (this.xc - this.xs == 0) {
				this.theta1 = 3.0 / 2.0 * Math.PI;
			} else {
				this.theta1 = Math.PI + Math.atan((this.yc - this.ys) / (this.xc - this.xs));
			}
		}
		if (this.xe < this.xc && this.ye < this.yc) {
			if (this.xc - this.xe == 0) {
				this.theta2 = 3.0 / 2.0 * Math.PI;
			} else {
				this.theta2 = Math.PI + Math.atan((this.yc - this.ye) / (this.xc - this.xe));
			} 
		}
		// Q4 -> 270°-360° 
		if (this.xs >= this.xc && this.ys < this.yc) {
			if (this.yc - this.ys == 0) {
				this.theta1 = 2.0 * Math.PI;
			} else {
				this.theta1 = 3.0 / 2.0 * Math.PI + Math.atan((this.xs - this.xc) / (this.yc - this.ys));
			}
		}
		if (this.xe >= this.xc && this.ye < this.yc) {
			if (this.yc - this.ye == 0) {
				this.theta2 = 2.0 * Math.PI;
			} else {
				this.theta2 = 3.0 / 2.0 * Math.PI + Math.atan((this.xe - this.xc) / (this.yc - this.ye));
			}
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
				double theta = theta1 + t1 * (theta2 - theta1);
				return yc + r * Math.sin(theta);
			}
		};
		
	}
	
	public static class Builder {
		
		private double xc;
		private double yc;
		private double xs;
		private double ys;
		private double xe;
		private double ye;

		public Builder xc(double xc) {
			this.xc = xc;
			return this;
		}

		public Builder yc(double yc) {
			this.yc = yc;
			return this;
		}

		public Builder xs(double xs) {
			this.xs = xs;
			return this;
		}

		public Builder ys(double ys) {
			this.ys = ys;
			return this;
		}

		public Builder xe(double xe) {
			this.xe = xe;
			return this;
		}

		public Builder ye(double ye) {
			this.ye = ye;
			return this;
		}

		public Arc build() {
			return new Arc(this);
		}
	}
	
}
