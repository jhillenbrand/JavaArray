package net.sytes.botg.array.geometry;

import net.sytes.botg.array.math.Mat;

/**
 * https://en.wikipedia.org/wiki/Gaussian_function
 * @author hillenbrand
 */
public class GaussianSurface extends SurfaceMesh {

	private double A = 1.0;
	private double x_0 = 0.0;
	private double y_0 = 0.0;
	private double sigma_x = 0.5;
	private double sigma_y = 0.5;
	
	public GaussianSurface() {
		
	}
	
	public GaussianSurface(double A, double x_0, double y_0, double sigma_x, double sigma_y) {
		this.A = A;
		this.x_0 = x_0;
		this.y_0 = y_0;
		this.sigma_x = sigma_x;
		this.sigma_y = sigma_y;
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
		// TODO Auto-generated method stub
		
	}
	
	private double gaussian2d(double x, double y) {
		return this.A * Math.exp(-1 * (Math.pow(x - this.x_0, 2) / 2 / Math.pow(this.sigma_x, 2) + Math.pow(y - this.y_0, 2) / 2 / Math.pow(this.sigma_y, 2)));
	}
}
