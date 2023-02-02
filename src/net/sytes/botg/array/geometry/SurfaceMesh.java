package net.sytes.botg.array.geometry;

public abstract class SurfaceMesh implements ISurfaceMesh {

	protected double[] x;
	protected double[] y;
	protected double[][] Z;
	
	public SurfaceMesh() {
	}
	
	public double[] getX() {
		return this.x;
	}

	public double[] getY() {
		return this.y;
	}

	public double[][] getZ() {
		return this.Z;
	}
	
}
