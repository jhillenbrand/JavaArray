package net.sytes.botg.array.geometry;

import net.sytes.botg.array.math.Vec;

public abstract class SurfaceMesh implements ISurfaceMesh {

	protected double[] x;
	protected double[] y;
	protected double[][] Z;
	
	public SurfaceMesh() {
	}
	
	@Override
	public void create() {
		double[] x = Vec.linspace(-1.0, 1.0, 100);
		double[] y = Vec.linspace(-1.0, 1.0, 100);
		this.create(x, y);
	}
	
	public double[] x() {
		return this.x;
	}

	public double[] y() {
		return this.y;
	}

	public double[][] z() {
		return this.Z;
	}
	
}
