package net.sytes.botg.array.math;

import net.sytes.botg.array.ArUtils;

public class Mat {
		
	// Suppress default constructor for noninstantiability
	private Mat() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
		
	/**
	 * computes the distance matrix for each point in points to all other points
	 * <br>here points[][] is a 2D double array, where the first dimension is the number of rows and the 2nd dimension the number of features/coordinates of the points
	 * <hr>Example 1:
	 * <br>if the distance between two Points P1 {1, 0} and P2 {1, 1} in 2D should be computed, points should look the following way:
	 * <br>point[0][0] = 1;
	 * <br>point[0][1] = 0;
	 * <br>point[1][0] = 1;
	 * <br>point[1][1] = 1;
	 * <hr>Example 2:
	 * <br>distance between two Points P1 {1, 1, 2} and P2 {0, 2, 3} in 3D
	 * <br>point[0][0] = 1;
	 * <br>point[0][1] = 1;
	 * <br>point[0][2] = 2;
	 * <br>point[1][0] = 0;
	 * <br>point[1][1] = 2;
	 * <br>point[1][2] = 3;
	 * <hr>The result returned is a symmetric matrix D_ij, that reads in the following way:
	 * <br>D_ij is the distance between Point i and Point j
	 * <br>all values D_ii are 0 (reads: distance from Point i to Point i --> 0)
	 * @param points
	 * @return
	 */
	public static double[][] distanceMatrix(double[][] points) {
		ArUtils.checkForNull(points);
		ArUtils.checkForEmpty(points);
		int n = points.length;
		double[][] D_ij = new double[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				// case 1 - skip diagonal
				if (i == j) {
					break;
				}
				
				// case 2 - symmetry
				if (j >= i) {
					break;
				}
				
				// case 3
				double d = Vec2Scalar.distance(points[i], points[j]);
				D_ij[i][j] = d;
				D_ij[j][i] = d;	// symmetry
				
				//ArUtils.print(D_ij);
				
			}
		}
		return D_ij;		
	}
	
	public static double[][] norm(double[][] x){
		int n = x.length;
		int m = x[0].length;		
		double xMin = Mat2Scalar.min(x);
		double xMax = Mat2Scalar.max(x);
		double[][] xn = new double[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				xn[i][j] = (x[i][j] - xMin) / (xMax - xMin);
			}
		}
		return xn;
	}
	
	/**
	 * normalization of vector {@code x} with zscore
	 * @param x
	 * @return
	 */
	public static double[][] zscore(double[][] x) {
		ArUtils.checkForEmpty(x);
		double mean = Mat2Scalar.mean(x);
		double sigma = Mat2Scalar.variance(x);
		int n = x.length;
		int m = x[0].length;
		double[][] z = new double[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				z[i][j] = (x[i][j] - mean) / sigma;
			}
		}
		return z;
	}	
}