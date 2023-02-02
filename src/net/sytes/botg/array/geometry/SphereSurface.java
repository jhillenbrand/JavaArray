package net.sytes.botg.array.geometry;

import net.sytes.botg.array.math.Mat;
import net.sytes.botg.array.math.Vec;

public class SphereSurface extends SurfaceMesh 	{

	private int n = 100;
	private double x_0 = 0.0;
	private double y_0 = 0.0;
	private double z_0 = 0.0;
	private double r = 1.0;
	
	public SphereSurface() {
		
	}
	
	public SphereSurface(int n, double r, double x_0, double y_0, double z_0) {
		this.n = n;
		this.r = r;
		this.x_0 = x_0;
		this.y_0 = y_0;
		this.z_0 = z_0;
	}
	
	/**
	 * 
	 */
	@Override
	public void create(double[] x, double[] y) {
		throw new UnsupportedOperationException(SphereSurface.class.getSimpleName() +  " is not defined for x[] and y[]");
	}

	@Override
	public void create() {
		double[] theta = Vec.linspace(-Math.PI, Math.PI, this.n);
		double[] phi = Vec.linspace(-Math.PI / 2, Math.PI / 2, this.n);
		
		double[] x = new double[this.n];
		double[] y = new double[this.n];
		double[][] Z = new double[this.n][this.n];
		
		for (int i = 0; i < this.n; i++) {
			x[i] = this.r * Math.sin(theta[i]) * Math.cos(phi[i]) + this.x_0;
			y[i] = this.r * Math.sin(theta[i]) * Math.cos(phi[i]) + this.y_0;			
		}
		
		double[][][] XY = Mat.meshgrid(x, y);
		double[][] X = XY[0];
		double[][] Y = XY[1];
		
		for (int i = 0; i < this.n; i++) {
			for (int j = 0; j < this.n; j++) {
				Z[i][j] = this.z_0 + Math.sqrt(Math.pow(this.r, 2) - Math.pow(X[i][j] - this.x_0, 2) - Math.pow(Y[i][j] - this.y_0, 2));
			}
		}
				
		this.x = x;
		this.y = y;
		
		this.Z = Z;
	}

}
