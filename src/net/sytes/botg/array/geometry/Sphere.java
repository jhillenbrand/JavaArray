package net.sytes.botg.array.geometry;

import net.sytes.botg.array.ConvertArray;
import net.sytes.botg.array.math.Vec;

/**
 * 3D Surface representin a semi sphere
 * <br>always insert two arrays t1 = [0 .. 1], t2 = [0 .. 1] with an arbitrary resolution in the create(...) method
 * @author hillenbrand
 *
 */
public class Sphere extends Surface3D {

	private double D;
	private double x_0;
	private double y_0;
	private double z_0;
	
	public Sphere() {
		this(new Builder());
	}
	
	private Sphere(Builder builder) {
		this.D = builder.D;
		this.x_0 = builder.x_0;
		this.y_0 = builder.y_0;
		this.z_0 = builder.z_0;
	}
	
	public static class Builder {
		private double D = 1.0;
		private double x_0 = 0.0;
		private double y_0 = 0.0;
		private double z_0 = 0.0;

		public Builder D(double D) {
			this.D = D;
			return this;
		}

		public Builder x_0(double x_0) {
			this.x_0 = x_0;
			return this;
		}

		public Builder y_0(double y_0) {
			this.y_0 = y_0;
			return this;
		}

		public Builder z_0(double z_0) {
			this.z_0 = z_0;
			return this;
		}

		public Sphere build() {
			return new Sphere(this);
		}
	}
	
	@Override
	public void create(double[] t1, double[] t2) {
		
		this.t1 = t1;
		this.t2 = t2;
		
		SemiSphere ss = new SemiSphere.Builder()
				.D(this.D)
				.x_0(this.x_0)
				.y_0(this.y_0)
				.z_0(this.z_0)
				.build();
		
		ss.create(t1, t2);
		
		this.x = ConvertArray.concat(ss.x(), ss.x());
		this.y = ConvertArray.concat(ss.y(), ss.y());
		this.z = ConvertArray.concat(ss.z(), Vec.product(ss.z(), -1.0));
		
	}
}
