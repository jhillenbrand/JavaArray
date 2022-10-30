package net.sytes.botg.array.geometry;

public abstract class Curve3D implements ICurve3D {
	
	protected double[] x;
	protected double[] y;
	protected double[] z;
	protected double[] t;
	
	protected FunctionHandle xHandle;
	protected FunctionHandle yHandle;
	protected FunctionHandle zHandle;
	
	public void create(double[] t) {
		this.t = t;
		this.x = new double[t.length];
		this.y = new double[t.length];
		this.z = new double[t.length];
		for (int i = 0; i < t.length; i++) {
			this.x[i] = this.xHandle.apply(this.t[i]);
			this.y[i] = this.yHandle.apply(this.t[i]);
			this.z[i] = this.zHandle.apply(this.t[i]);
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
	
	
	public double[] t() {
		return this.t;
	}
}
