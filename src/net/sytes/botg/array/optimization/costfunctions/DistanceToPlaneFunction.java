package net.sytes.botg.array.optimization.costfunctions;

import net.sytes.botg.array.optimization.ObjectFunction;

public class DistanceToPlaneFunction extends ObjectFunction {

	public DistanceToPlaneFunction() {
		super(4);
	}

	@Override
	public double apply(double[][] X_ij) {
		double sum_squared = 0.0;
		for (int i = 0; i < X_ij.length; i++) {
			double d = this.apply(X_ij[i]);
			sum_squared = sum_squared + d * d;
		}
		return sum_squared;
	}

	/**
	 * equation taken from 
	 * <a href="https://mathinsight.org/distance_point_plane">Link</a>
	 */
	@Override
	public double apply(double[] x_j) {
		if (x_j.length != 3) {
			throw new IllegalArgumentException("Vector x_j must be of size 3");
		}
		double A = this.theta_j[0];
		double B = this.theta_j[1];
		double C = this.theta_j[2]; 
		double D = this.theta_j[3];
		
		double d = Math.abs(A * x_j[0] + B * x_j[1] + C * x_j[2] + D) / Math.sqrt(A * A + B * B + C * C);
		
		return d;
	}

}
