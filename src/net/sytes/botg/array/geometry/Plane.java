package net.sytes.botg.array.geometry;

public class Plane extends Surface3D {

	private double nx;
	private double ny;
	private double nz;
	private double x0;
	private double y0;
	private double z0;
	private double xmin;
	private double xmax;
	private double ymin;
	private double ymax;

	private Plane(Builder builder) {
		this.nx = builder.nx;
		this.ny = builder.ny;
		this.nz = builder.nz;
		this.x0 = builder.x0;
		this.y0 = builder.y0;
		this.z0 = builder.z0;
		this.xmin = builder.xmin;
		this.xmax = builder.xmax;
		this.ymin = builder.ymin;
		this.ymax = builder.ymax;
		
		// create function handles
		this.xHandle = new FunctionHandle() {
			@Override
			public double apply(double t1, double t2) {
				return xmin + t1 * (xmax - xmin);
			}
		};
				
		this.yHandle = new FunctionHandle() {
			@Override
			public double apply(double t1, double t2) {
				return ymin + t2 * (ymax - ymin);
			}
		};
				
		this.zHandle = new FunctionHandle() {
			@Override
			public double apply(double t1, double t2) {
				double x = xmin + t1 * (xmax - xmin);
				double y = ymin + t2 * (ymax - ymin);
				double z = z0 - (nx * (x - x0) + ny * (y - y0)) / nz;
				return z;
			}
		};
		
	}
	
	public static class Builder {
		private double nx;
		private double ny;
		private double nz;
		private double x0;
		private double y0;
		private double z0;
		private double xmin;
		private double xmax;
		private double ymin;
		private double ymax;

		public Builder nx(double nx) {
			this.nx = nx;
			return this;
		}

		public Builder ny(double ny) {
			this.ny = ny;
			return this;
		}

		public Builder nz(double nz) {
			this.nz = nz;
			return this;
		}

		public Builder x0(double x0) {
			this.x0 = x0;
			return this;
		}

		public Builder y0(double y0) {
			this.y0 = y0;
			return this;
		}

		public Builder z0(double z0) {
			this.z0 = z0;
			return this;
		}

		public Builder xmin(double xmin) {
			this.xmin = xmin;
			return this;
		}

		public Builder xmax(double xmax) {
			this.xmax = xmax;
			return this;
		}

		public Builder ymin(double ymin) {
			this.ymin = ymin;
			return this;
		}

		public Builder ymax(double ymax) {
			this.ymax = ymax;
			return this;
		}

		public Plane build() {
			return new Plane(this);
		}
	}
}
