package net.sytes.botg.array.optimization.algorithms;

import net.sytes.botg.array.optimization.Algorithm;
import net.sytes.botg.array.optimization.ObjectFunction;

public class MultipleLinearRegression extends Algorithm {

	/**
	 * learn rate
	 */
	private double alpha = 0.1;
	
	/**
	 * number of features
	 */
	private int f = 2;
	
	/**
	 * weights
	 */
	private double[] w;
	
	/**
	 * bias
	 */
	private double b;
	
	public MultipleLinearRegression(ObjectFunction objectFunction, ObjectFunction gradientFunction) {
		super(objectFunction, gradientFunction);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void run(double[][] X_ij) {
		// TODO Auto-generated method stub
	}
	
	
}
