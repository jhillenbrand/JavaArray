package net.sytes.botg.array.optimization;

public abstract class ObjectFunction implements IObjectFunction {
	
	protected int n;	
	protected double[] theta_j;
	
	public ObjectFunction(int numOfWeights) {
		this.n = numOfWeights;
		this.theta_j = new double[this.n];
	}
	
}
