package net.sytes.botg.array.geometry;

public class Helix extends Curve3D {

	private double D;
	private double P;
	private double z_0;
	private boolean right;

	public Helix() {
		this(new Builder());
	}
	
	private Helix(Builder builder) {
		this.D = builder.D;
		this.P = builder.P;
		this.z_0 = builder.z_0;
		this.right = builder.right;
		if (right) {
			this.xHandle = new FunctionHandle() {
				@Override
				public double apply(double t) {
					return D / 2.0 * Math.cos(2 * Math.PI * t);
				}
			};

			this.yHandle = new FunctionHandle() {
				@Override
				public double apply(double t) {
					return D / 2.0 * Math.sin(2 * Math.PI * t);
				}
			};

			this.zHandle = new FunctionHandle() {
				@Override
				public double apply(double t) {
					return P * t + z_0;
				}
			};
		} else {
			this.xHandle = new FunctionHandle() {
				@Override
				public double apply(double t) {
					return D / 2.0 * Math.cos(2 * Math.PI * t);
				}
			};

			this.yHandle = new FunctionHandle() {
				@Override
				public double apply(double t) {
					return - D / 2.0 * Math.sin(2 * Math.PI * t);
				}
			};

			this.zHandle = new FunctionHandle() {
				@Override
				public double apply(double t) {
					return P * t + z_0;
				}
			};
		}
	}
	
	public static class Builder {
		
		private double D = 10.0;
		private double P = 1.0;
		private double z_0 = 0.0;
		private boolean right = true;

		public Builder D(double D) {
			this.D = D;
			return this;
		}

		public Builder P(double P) {
			this.P = P;
			return this;
		}

		public Builder z_0(double z_0) {
			this.z_0 = z_0;
			return this;
		}

		public Builder right(boolean right) {
			this.right = right;
			return this;
		}

		public Helix build() {
			return new Helix(this);
		}
	}
}
