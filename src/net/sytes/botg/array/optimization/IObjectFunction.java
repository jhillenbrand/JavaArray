package net.sytes.botg.array.optimization;

public interface IObjectFunction {

	public double apply(double[] x_j);
	
	public double apply(double[][] X_ij);
		
}
