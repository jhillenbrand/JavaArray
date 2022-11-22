package net.sytes.botg.array.math;

import net.sytes.botg.array.ArUtils;

/**
 * Class that contains matrix to matrix operations
 * <br>if not specified in the method description, the matrices adhere the following dimension nomenclature:
 * <br>>> double[][] X = new double[m][n],
 * <br>where n is the number of rows and m the number of columns
 * @author hillenbrand
 *
 */
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
	 * <br>all values D_ii are 0 (reads: distance from Point i to Point i --&gt; 0)
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
	
	/**
	 * assuming {@code XYZ} is a matrix with 3 vectors of x-,y-,z-values of a position vector,
	 * <br>this methods returns the normed vectors
	 * @param XYZ
	 * @return
	 */
	public static double[][] normVectors(double[][] XYZ){
		ArUtils.checkForNull(XYZ);
		ArUtils.checkForEmpty(XYZ);
		int n = XYZ.length;
		int m = XYZ[0].length;
		double[][] XYZ_n = new double[n][m];
		for (int j = 0; j < m; j++) {
			double squaredSum = 0.0;
			for (int i = 0; i < n; i++) {
				squaredSum = squaredSum + XYZ[i][j] * XYZ[i][j]; 
			}
			double vecLen = Math.sqrt(squaredSum);
			for (int i = 0; i < n; i++) {
				XYZ_n[i][j] = XYZ[i][j] / vecLen;
			}
		}
		return XYZ_n;
	}
	
	/**
	 * normalize matrix entries in range from 0 to 1
	 * @param X
	 * @return
	 */
	public static double[][] normalize(double[][] X){
		int n = X.length;
		int m = X[0].length;		
		double xMin = Mat2Scalar.min(X);
		double xMax = Mat2Scalar.max(X);
		double[][] xn = new double[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				xn[i][j] = (X[i][j] - xMin) / (xMax - xMin);
			}
		}
		return xn;
	} 
	
	/**
	 * normalization of vector {@code x} with zscore
	 * @param X
	 * @return
	 */
	public static double[][] zscore(double[][] X) {
		ArUtils.checkForEmpty(X);
		double mean = Mat2Scalar.mean(X);
		double sigma = Mat2Scalar.variance(X);
		int n = X.length;
		int m = X[0].length;
		double[][] z = new double[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				z[i][j] = (X[i][j] - mean) / sigma;
			}
		}
		return z;
	}
	
	public static double[][] prod(double[][] X, double[][] Y){
		ArUtils.checkForNull(X, Y);
		checkMatrixProdDimensions(X, Y);
		
		int n1 = X.length;
		int m1 = X[0].length;
		
		
		return null;
	}
	
	/**
	 * combine two matrices based on their dimensions
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[][] appendByColumn(double[][] X, double[][] Y){
		int n1 = X[0].length;
		int m1 = X.length;
		int n2 = Y[0].length;
		int m2 = Y.length;
		if (n1 == n2) {
			double[][] Z = new double[m1 + m2][n1];
			int c;
			for (c = 0; c < m1; c++) {
				Z[c] = X[c];
			}
			int c2 = c;
			for (c = 0; c < m2; c++) {
				Z[c2] = Y[c];
				c2 = c2 + c;
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + ArUtils.matrixDimensionsToString(X) + " and Y" + ArUtils.matrixDimensionsToString(Y) + " do not have matching rows for appending by column!");
		}
	}
	
	private static  void checkMatrixProdDimensions(double[][] X, double[][] Y) {
		if (X[0].length != Y.length) {
			throw new IllegalArgumentException("Dimension mismatch for matrix multiplication");
		}
	}
	
	/**
	 * check if both matrices have same dimensions
	 * @param X
	 * @param Y
	 */
	public static void checkForMatchingDimensions(double[][] X, double[][] Y) {
		if (X.length != Y.length || X[0].length != Y[0].length) {
			throw new IllegalArgumentException("Dimension mismatch of matrices");
		}
	}
}
