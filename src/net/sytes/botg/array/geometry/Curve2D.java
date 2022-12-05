package net.sytes.botg.array.geometry;

public abstract class Curve2D extends Curve3D {
	
	@Override
	public void create(double[] t) {
		this.t = t;
		this.x = new double[t.length];
		this.y = new double[t.length];
		for (int i = 0; i < t.length; i++) {
			this.x[i] = this.xHandle.apply(this.t[i]);
			this.y[i] = this.yHandle.apply(this.t[i]);
		}	
	}
	
	/**
	 * not applicable for {@code Curve2D}
	 */
	@Override
	public double[] z() {
		throw new UnsupportedOperationException("Method z() cannot be executed for " + Curve2D.class.getSimpleName() + " objects");
	}
}
