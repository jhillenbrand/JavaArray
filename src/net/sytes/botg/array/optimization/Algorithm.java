package net.sytes.botg.array.optimization;

public abstract class Algorithm implements IAlgorithm {

	protected ObjectFunction objectFunction;
	protected ObjectFunction gradientFunction;
	
	public Algorithm(ObjectFunction objectFunction) {
		this(objectFunction, null);
	}
	
	public Algorithm(ObjectFunction objectFunction, ObjectFunction gradientFunction) {
		this.objectFunction = objectFunction;
		this.gradientFunction = gradientFunction;
	}
	
	public ObjectFunction createObjectFunction(double[][] X_ij) {
		
		this.run(X_ij);		
		return this.objectFunction;		
	}
	
	private void verifyInputs(double[][] X_ij) {
		 
	}
	
}
