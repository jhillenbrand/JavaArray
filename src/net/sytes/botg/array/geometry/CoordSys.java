package net.sytes.botg.array.geometry;

import net.sytes.botg.array.math.Mat;
import net.sytes.botg.array.math.Vec;

/**
 * Class representing a coordinate system defined by its base vectors ({@code base} in a orthogonal euclidean x, y, z system
 * and its {@code origin} vector
 * <br>
 * <img src="CoordSys.png"/>
 * @author hillenbrand
 *
 */
public class CoordSys {

	private double[][] base = new double[][] {{1, 0, 0},{0, 1, 0},{0, 0, 1}};
	private double[] origin = new double[] {0, 0, 0};	
	private String[] labels = new String[] {"x", "y", "z"};
	
	public static CoordSys defaultSystem() {
		return new CoordSys();
	}
	
	public void origin(double[] o) {
		if (o.length < 2 || o.length > 3) {	
			throw new IllegalArgumentException(CoordSys.class.getSimpleName() + " is only defined for 2D or 3D coordinates, specify the origin accordingly.");
		}
		this.origin = o;
	}
	
	public double[] origin() {
		return this.origin;
	}
	
	/**
	 * moves the {@code origin} by vector {@code b}
	 * @param b
	 */
	public void move(double[] b) {
		if (b.length < 2 || b.length > 3) {	
			throw new IllegalArgumentException(CoordSys.class.getSimpleName() + " is only defined for 2D or 3D coordinates, specify the translation vector b accordingly.");
		}
		this.origin = Vec.plus(this.origin, b);
	}
	
	/**
	 * sets the new {@code base}
	 * @param base
	 */
	public void base(double[][] base) {
		if (base.length < 2 || base.length > 3) {	
			throw new IllegalArgumentException(CoordSys.class.getSimpleName() + " is only defined for 2D or 3D coordinates, specify the base matrix accordingly.");
		}
		this.base = base;
	}
	
	public double[][] base(){
		return this.base;
	}
	
	public void labels(String[] labels) {
		if (labels.length < 2 || labels.length > 3) {	
			throw new IllegalArgumentException(CoordSys.class.getSimpleName() + " is only defined for 2D or 3D coordinates, specify the labels accordingly.");
		}
		this.labels = labels;
	}
	
	public String[] labels() {
		return this.labels;
	}
	
	/**
	 * rotates the {@code base} coordinate system by rotation matrix {@code R}
	 * @param R
	 */
	public void rot(double[][] R) {
		if (R.length < 2 || R.length > 3) {	
			throw new IllegalArgumentException(CoordSys.class.getSimpleName() + " is only defined for 2D or 3D coordinates, specify the rotation matrix R accordingly.");
		}
		this.base = Mat.product(this.base, R);
	}
	
	/**
	 * rotates the {@code base} around 1st axis by angle {@code phi}
	 * @param phi
	 */
	public void rot1(double phi) {
		double[][] R = new double[][] {
			{ 1,                              0,                             0},
			{ 0,  Math.cos(Math.toRadians(phi)), Math.sin(Math.toRadians(phi))},
			{ 0, -Math.sin(Math.toRadians(phi)), Math.cos(Math.toRadians(phi))}
		};
	
		this.base = Mat.product(this.base, R);
	}
	
	/**
	 * rotates the {@code base} around 2nd axis by angle {@code phi}
	 * @param phi
	 */
	public void rot2(double phi) {
		double[][] R = new double[][] {
			{ Math.cos(Math.toRadians(phi)),                              0, -Math.sin(Math.toRadians(phi))},
			{                             0,                              1,                              0},
			{ Math.sin(Math.toRadians(phi)),                              0,  Math.cos(Math.toRadians(phi))}
		};
	
		this.base = Mat.product(this.base, R);
	}
	
	/**
	 * rotates the {@code base} around 3rd axis by angle {@code phi}
	 * @param phi
	 */
	public void rot3(double phi) {
		double[][] R = new double[][] {
				{ Math.cos(Math.toRadians(phi)), Math.sin(Math.toRadians(phi)), 0},
				{-Math.sin(Math.toRadians(phi)), Math.cos(Math.toRadians(phi)), 0},
				{                             0,                             0, 1}
			};
		
		this.base = Mat.product(this.base, R);		
	}
	
	/**
	 * transforms a vector {@code x} in the euclidean base system {{1,0,0}, {0,1,0}, {0,0,1}} (or a previous transformed system) into a vector {@code x'} of this coordinate system
	 * @param x
	 * @return
	 */
	public double[] transform(double[] x) {
		double[] x_b = Vec.minus(x, this.origin);
		double[][] X_B = Vec.matrix(x_b, 1);
		double[][] X_ = Mat.product(this.base, X_B);
		return X_[0];
	}
	
	@Override
	public CoordSys clone() {
		CoordSys cos = new CoordSys();
		cos.base = Mat.copy(this.base);
		cos.origin = Vec.copy(this.origin);
		return cos;
	}
	
}
