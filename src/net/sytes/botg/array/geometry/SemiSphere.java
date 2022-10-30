package net.sytes.botg.array.geometry;

/**
 * 3D Surface representin a semi sphere
 * <br>always insert two arrays t1 = [0 .. 1], t2 = [0 .. 1] with an arbitrary resolution in the create(...) method
 * @author hillenbrand
 *
 */
public class SemiSphere extends Surface3D {

	private double D;
	private double x_0;
	private double y_0;
	private double z_0;
	private boolean top;

	public SemiSphere() {
		this(new Builder());
	}

	private SemiSphere(Builder builder) {
		this.D = builder.D;
		this.x_0 = builder.x_0;
		this.y_0 = builder.y_0;
		this.z_0 = builder.z_0;
		this.top = builder.top;
		
		// create function handles
		this.xHandle = new FunctionHandle() {
			@Override
			public double apply(double t1, double t2) {
				double theta = t1 * 2 * Math.PI;
				double phi = t2 * Math.PI / 2;
				return x_0 + D / 2.0 * Math.cos(theta) * Math.sin(phi);
			}
		};
		
		this.yHandle = new FunctionHandle() {
			@Override
			public double apply(double t1, double t2) {
				double theta = t1 * 2 * Math.PI;
				double phi = t2 * Math.PI / 2;
				return y_0 + D / 2.0 * Math.sin(theta) * Math.sin(phi);
			}
		};
		
		this.zHandle = new FunctionHandle() {
			@Override
			public double apply(double t1, double t2) {
				double theta = t1 * 2 * Math.PI;
				double phi = t2 * Math.PI / 2;
				if (top) {
					return z_0 + D / 2.0 * Math.cos(phi);
				} else {
					return z_0 - D / 2.0 * Math.cos(phi);
				}
			}
		};
	}
	
	public static class Builder {
		private double D = 1.0;
		private double x_0 = 0.0;
		private double y_0 = 0.0;
		private double z_0 = 0.0;
		private boolean top = true;

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

		public Builder z_0(double z_0) {
			this.z_0 = z_0;
			return this;
		}
		
		public Builder top(boolean top) {
			this.top = top;
			return this;
		}

		public SemiSphere build() {
			return new SemiSphere(this);
		}
	}
}
