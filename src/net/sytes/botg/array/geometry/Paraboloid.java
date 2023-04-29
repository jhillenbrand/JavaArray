package net.sytes.botg.array.geometry;

import net.sytes.botg.array.math.Mat;
import net.sytes.botg.array.math.Vec;

public class Paraboloid extends SurfaceMesh {

	private double a;
	private double b;
	private double z0;
	private boolean elliptic;

	public Paraboloid() {
		this(new Builder());
	}

	private Paraboloid(Builder builder) {		
		this.a = builder.a;
		this.b = builder.b;
		this.z0 = builder.z0;
		this.elliptic = builder.elliptic;
	}

	public static class Builder {
		
		private double a = 1.0;
		private double b = 1.0;
		private double z0 = 0.0;
		private boolean elliptic = true;
		
		public Builder a(double a) {
			this.a = a;
			return this;
		}

		public Builder b(double b) {
			this.b = b;
			return this;
		}

		public Builder z0(double z0) {
			this.z0 = z0;
			return this;
		}

		public Builder elliptic(boolean elliptic) {
			this.elliptic = elliptic;
			return this;
		}

		public Paraboloid build() {
			return new Paraboloid(this);
		}
	}

	@Override
	public void create(double[] x, double[] y) {
		int n = x.length;
		int m = y.length;
		this.x = x;
		this.y = y;
		this.Z = new double[n][m];
		double[][][] XY = Mat.meshgrid(x, y);
		double[][] X = XY[0];
		double[][] Y = XY[1];
		if (this.elliptic) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					this.Z[i][j] = this.paraboloid(X[i][j], Y[i][j]);
				}
			}
		} else {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					this.Z[i][j] = this.hyperboloid(X[i][j], Y[i][j]);
				}
			}
		}
	}

	@Override
	public double create(double x, double y) {
		if (this.elliptic) {
			return this.paraboloid(x, y);
		} else {
			return this.hyperboloid(x, y);
		}
	}

	@Override
	public void create() {
		double[] x = Vec.linspace(-1.0, 1.0, 100);
		double[] y = Vec.linspace(-1.0, 1.0, 100);
		this.create(x, y);
	}
	
	private double hyperboloid(double x, double y) {
		return this.z0 + y * y / this.b / this.b - x * x / this.a / this.a;
	}
	
	private double paraboloid(double x, double y) {
		return this.z0 + y * y / this.b / this.b + x * x / this.a / this.a;
	}
	
}
