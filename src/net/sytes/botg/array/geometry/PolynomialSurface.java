package net.sytes.botg.array.geometry;

import net.sytes.botg.array.math.Mat;

/**
 * defines a polynomial surface
 * <br>f(x,y) = c<sub>0</sub> + c<sub>1</sub>&middot;x + c<sub>1</sub>&middot;x + c<sub>2</sub>&middot;y + c<sub>3</sub>&middot;x&middot;y + c<sub>4</sub>&middot;x<sup>2</sup> + c<sub>5</sub>&middot;y<sup>2</sup> + c<sub>6</sub>&middot;x<sup>2</sup>&middot;y + c<sub>7</sub>&middot;x&middot;y<sup>2</sup> + c<sub>8</sub>&middot;x<sup>3</sup> + c<sub>9</sub>&middot;y<sup>3</sup> + ...
 * <br>where c<sub>i</sub> are the polynomial weights
 * <br>or defined as a matrix with coefficients a<sub>ij</sub>:
 * <br>f(x,y) = &sum;<sub>i=0</sub>&sum;<sub>j=0</sub>a<sub>ij</sub>&middot;x<sup>i</sup>&middot;y<sup>j</sup>
 * @author hillenbrand
 */
public class PolynomialSurface extends SurfaceMesh {

	private double[][] a_ij;

	public PolynomialSurface() {
		this(new Builder());
	}

	private PolynomialSurface(Builder builder) {
		this.a_ij = builder.a_ij;
	}

	public static class Builder {
		private double[][] a_ij;

		public Builder a_ij(double[][] a_ij) {
			this.a_ij = a_ij;
			return this;
		}

		public Builder weights(double[] weights) {
			int n = weights.length;
			int s = (int) Math.ceil(Math.sqrt(n));
			this.a_ij = new double[s][s];
			int k = 0;
			int i = 0;
			int j = 0;
			while (k < n) {
				this.a_ij[i][j] = weights[k];
				++i;
				++k;
				if (k < n) {
					this.a_ij[i][j] = weights[k];
				}
				--i;
				++j;
				++k;
				if (k < n) {
					this.a_ij[i][j] = weights[k];
				}
				++i;
				++k;
			}
			return this;
		}
		
		public PolynomialSurface build() {
			return new PolynomialSurface(this);
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
				this.Z[i][j] = this.polySurface(X[i][j], Y[i][j]);
			}
		}
	}

	private double polySurface(double x, double y) {
		double p = 0.0;

		// assemble x-/y-matrices
		int n = this.a_ij.length;
		double[][] X = new double[1][n];
		double[][] Y = new double[n][1];

		for (int i = 0; i < n; i++) {
			X[0][i] = Math.pow(x, i);
			Y[i][0] = Math.pow(y, i);
		}
		
		// matrix multiplication
		
		double[][] b_ij = Mat.product(this.a_ij, Y);
		double[][] c_ij = Mat.product(X, b_ij);
		
		return c_ij[0][0];
	}
}
