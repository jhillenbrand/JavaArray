package net.sytes.botg.array.geometry;

import net.sytes.botg.array.math.Mat;
import net.sytes.botg.array.math.Vec;

/**
 * f(x, y) = A &middot; exp(-((x-x<sub>0</sub>)<sup>2</sup>/(2&middot;&sigma;<sub>x</sub><sup>2</sup>)+(y-y<sub>0</sub>)<sup>2</sup>/(2&middot;&sigma;<sub>y</sub><sup>2</sup>)))
 * <br>https://en.wikipedia.org/wiki/Gaussian_function
 * @author hillenbrand
 */
public class GaussianSurface extends SurfaceMesh {

	private double A;
	private double x_0;
	private double y_0;
	private double sigma_x;
	private double sigma_y;

	public GaussianSurface() {
		this(new Builder());
	}

	private GaussianSurface(Builder builder) {
		this.A = builder.A;
		this.x_0 = builder.x_0;
		this.y_0 = builder.y_0;
		this.sigma_x = builder.sigma_x;
		this.sigma_y = builder.sigma_y;
	}

	public static class Builder {
		private double A = 1.0;
		private double x_0 = 0.0;
		private double y_0 = 0.0;
		private double sigma_x = 0.5;
		private double sigma_y = 0.5;

		public Builder A(double A) {
			this.A = A;
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

		public Builder sigma_x(double sigma_x) {
			this.sigma_x = sigma_x;
			return this;
		}

		public Builder sigma_y(double sigma_y) {
			this.sigma_y = sigma_y;
			return this;
		}

		public GaussianSurface build() {
			return new GaussianSurface(this);
		}
	}

	@Override
	public void create(double[] x, double[] y) {
		int n = x.length;
		int m = y.length;
		this.Z = new double[n][m];
		double[][][] XY = Mat.meshgrid(x, y);
		double[][] X = XY[0];
		double[][] Y = XY[1];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				this.Z[i][j] = this.gaussian2d(X[i][j], Y[i][j]);
			}
		}
	}

	@Override
	public void create() {
		double[] x = Vec.linspace(-1.0, 1.0, 100);
		double[] y = Vec.linspace(-1.0, 1.0, 100);
		this.create(x, y);
	}

	private double gaussian2d(double x, double y) {
		return this.A * Math.exp(-1 * (Math.pow(x - this.x_0, 2) / 2 / Math.pow(this.sigma_x, 2) + Math.pow(y - this.y_0, 2) / 2 / Math.pow(this.sigma_y, 2)));
	}

	@Override
	public double create(double x, double y) {
		return this.gaussian2d(x, y);
	}
	
	public double[] weights() {
		return new double[] { this.A, this.x_0, this.y_0, this.sigma_x, this.sigma_y};
	}
}
