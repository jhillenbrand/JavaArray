package net.sytes.botg.array.geometry;

public abstract class Surface3D implements ISurface3D {

	protected double[] x;
	protected double[] y;
	protected double[] z;
	protected double[] t1;
	protected double[] t2;
	
	protected FunctionHandle xHandle;
	protected FunctionHandle yHandle;
	protected FunctionHandle zHandle;
		
	public void create(double[] t1, double[] t2) {
		this.t1 = t1;
		this.t2 = t2;
		this.x = new double[t1.length * t2.length];
		this.y = new double[t1.length * t2.length];
		this.z = new double[t1.length * t2.length];
		int c = 0;
		for (int i = 0; i < t1.length; i++) {
			for (int j = 0; j < t2.length; j++) {				
				this.x[c] = this.xHandle.apply(this.t1[i], this.t2[j]);
				this.y[c] = this.yHandle.apply(this.t1[i], this.t2[j]);
				this.z[c] = this.zHandle.apply(this.t1[i], this.t2[j]);
				++c;
			}
		}	
	}
	
	public double[] x() {
		return this.x;
	}
	
	public double[] y() {
		return this.y;
	}
	
	public double[] z() {
		return this.z;
	}
		
	public double[] t1() {
		return this.t1;
	}
	
	public double[] t2() {
		return this.t2;
	}
	
}
