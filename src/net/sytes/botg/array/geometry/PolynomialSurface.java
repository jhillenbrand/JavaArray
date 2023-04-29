package net.sytes.botg.array.geometry;

import net.sytes.botg.array.math.Mat;

/**
 * defines a polynomial surface
 * <br>f(x,y) = c<sub>0</sub> + c<sub>1</sub>&middot;x + c<sub>1</sub>&middot;x + c<sub>2</sub>&middot;y + c<sub>3</sub>&middot;x&middot;y + c<sub>4</sub>&middot;x<sup>2</sup> + c<sub>5</sub>&middot;y<sup>2</sup> + c<sub>6</sub>&middot;x<sup>2</sup>&middot;y + c<sub>7</sub>&middot;x&middot;y<sup>2</sup> + c<sub>8</sub>&middot;x<sup>3</sup> + c<sub>9</sub>&middot;y<sup>3</sup> + ...
 * <br>where c<sub>i</sub> are the polynomial weights
 * <br>or defined as a matrix with coefficients a<sub>ij</sub>:
 * <br>f(x,y) = &sum;<sub>i=0</sub>&sum;<sub>j=0</sub>a<sub>ij</sub>&middot;x<sup>i</sup>&middot;y<sup>j</sup>
 * @author hillenbrand
 */
public class PolynomialSurface extends SurfaceMesh {

	private double[][] a_ij;

	public PolynomialSurface() {
		this(new Builder());
	}

	private PolynomialSurface(Builder builder) {
		this.a_ij = builder.a_ij;
	}

	public static class Builder {
		private double[][] a_ij;

		public Builder a_ij(double[][] a_ij) {
			this.a_ij = a_ij;
			return this;
		}

		public Builder weights(double[] weights) {
			this.a_ij = matrix(weights);
			return this;
		}
		
		public PolynomialSurface build() {
			return new PolynomialSurface(this);
		}
	}

	@Override
	public void create(double[] x, double[] y) {
		int n = x.length;
		int m = y.length;
		this.Z = new double[n][m];
		double[][][] XY = Mat.meshgrid(x, y);
		double[][] X = XY[0];
		double[][] Y = XY[1];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				this.Z[i][j] = this.polySurface(X[i][j], Y[i][j]);
			}
		}
	}

	@Override
	public double create(double x, double y) {		
		return this.polySurface(x, y);
	}

	private double polySurface(double x, double y) {
		
		// assemble x-/y-matrices
		int n = this.a_ij.length;
		double[][] X = new double[1][n];
		double[][] Y = new double[n][1];

		for (int i = 0; i < n; i++) {
			X[0][i] = Math.pow(x, i);
			Y[i][0] = Math.pow(y, i);
		}
		
		// matrix multiplication
		
		double[][] b_ij = Mat.product(this.a_ij, Y);
		double[][] c_ij = Mat.product(X, b_ij);
		
		return c_ij[0][0];
	}
	
	private static double[][] matrix(double[] weights) {
		int n = weights.length;
		int s = (int) Math.ceil(Math.sqrt(n));
		double[][] a_ij = new double[s][s];
		int k = 0;
		int i = 0;
		int j = 0;
		int p = a_ij.length - 1;
		int nd = 2 * p + 1;
		boolean stop = false;
		
		// go through diagonal matrix elements
		for (int d = 0; d < nd; d++) {
			
			// even diagonals
			if (d % 2 == 0) {
				
				a_ij[i][j] = weights[k];
				++k;
				if (k >= n) {
					break;
				}
				// go through elements on diagonal
				int ii = i;
				int jj = j;
				while (ii + 1 <= p && jj - 1 >= 0) {					
					ii = ii + 1;
					jj = jj - 1;
					a_ij[ii][jj] = weights[k];
					++k;
					if (k >= n) {
						stop = true;
						break;
					}
					// set transpose element
					a_ij[jj][ii] = weights[k];
					++k;
					if (k >= n) {
						stop = true;
						break;
					}
				}
				
				++i;
				
			// odd diagonals
			} else {
				
				a_ij[i][j] = weights[k];
				++k;
				if (k >= n) {
					break;
				}
				// set transpose element
				a_ij[j][i] = weights[k];
				++k;
				if (k >= n) {
					break;
				}
				// go through elements on diagonal
				int ii = i;
				int jj = j;
				while (ii + 1 <= p && jj - 1 >= 0) {					
					ii = ii + 1;
					jj = jj - 1;
					a_ij[ii][jj] = weights[k];
					++k;
					if (k >= n) {
						stop = true;
						break;
					}
					// set transpose element
					a_ij[jj][ii] = weights[k];
					++k;
					if (k >= n) {
						stop = true;
						break;
					}
				}
				
				++j;
			}
			
			if (stop) {
				break;
			}			
		}
		return a_ij;
	}
	
	public double[] weights() {		
		int n = Mat.numel(this.a_ij);
		int k = 0;
		int i = 0;
		int j = 0;
		int p = this.a_ij.length - 1;
		int nd = 2 * p + 1;
		boolean stop = false;
		double[] weights = new double[n];
		
		// go through diagonal matrix elements
		for (int d = 0; d < nd; d++) {
			
			// even diagonals
			if (d % 2 == 0) {
				
				weights[k] = this.a_ij[i][j];
				++k;
				if (k >= n) {
					break;
				}
				// go through elements on diagonal
				int ii = i;
				int jj = j;
				while (ii + 1 <= p && jj - 1 >= 0) {					
					ii = ii + 1;
					jj = jj - 1;
					weights[k] = this.a_ij[ii][jj];
					++k;
					if (k >= n) {
						stop = true;
						break;
					}
					// set transpose element
					weights[k] = this.a_ij[jj][ii];
					++k;
					if (k >= n) {
						stop = true;
						break;
					}
				}
				
				++i;
				
			// odd diagonals
			} else {
				
				weights[k] = this.a_ij[i][j];
				++k;
				if (k >= n) {
					break;
				}
				// set transpose element
				weights[k] = this.a_ij[j][i];
				++k;
				if (k >= n) {
					break;
				}
				// go through elements on diagonal
				int ii = i;
				int jj = j;
				while (ii + 1 <= p && jj - 1 >= 0) {					
					ii = ii + 1;
					jj = jj - 1;
					weights[k] = this.a_ij[ii][jj];
					++k;
					if (k >= n) {
						stop = true;
						break;
					}
					// set transpose element
					weights[k] = this.a_ij[jj][ii];
					++k;
					if (k >= n) {
						stop = true;
						break;
					}
				}
				
				++j;
			}
			
			if (stop) {
				break;
			}			
		}

		
		/*
		while (k < n) {		
			// set weight for current depth
			weights[k] = this.a_ij[i + d][j];
			++k;
			
			// go one deeper
			++d;
			if (k < n) {
				weights[k] = this.a_ij[i + d][j];
			}			
			++k;
			// set weight for transposed element
			if (k < n) {
				weights[k] = this.a_ij[i][j + d];
			}			
			++k;
		}
		*/		
		return weights;
	}
	
	public double[][] a_ij(){
		return this.a_ij;
	}
	
}
