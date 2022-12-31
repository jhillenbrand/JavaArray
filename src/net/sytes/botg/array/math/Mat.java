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
	
	/**
	 * returns the matrix multiplication of {@code X} and {@code Y}
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[][] prod(double[][] X, double[][] Y){
		ArUtils.checkForNull(X, Y);
		checkMatrixProdDimensions(X, Y);
		
        double[][] Z = new double[X.length][Y[0].length];
        for (int i = 0; i < X.length; i++) {
            for (int j = 0; j < Y[0].length; j++) {
                double sum = 0;
                for (int k = 0; k < X[i].length; k++) {
                    sum += X[i][k] * Y[k][j];
                }
                Z[i][j] = sum;
            }
        }
        return Z;
	}
	
	/**
	 * solve the linear equation system y<sub>j</sub> = A<sub>ij</sub> x<sub>i</sub> for x<sub>i</sub>
	 * @param A
	 * @param y
	 * @return
	 */
	public static double[] linSolve(double[][] A, double[] y) {
		
		if (A.length != y.length) {
			throw new IllegalArgumentException("First matrix dimension of A (" + A.length + ") and length of vector y (" + y.length + ") must be equal");
		}
		
		int n = A.length;
		
		// put y-vector into matrix Y for matrix multiplication
		double[][] Y = Vec2Mat.vec2Mat(y, n, 1, false);
		
		double[][] X = Mat.prod(Mat.inverse(A), Y); 
		
		// reshape the result matrix back to vector
		double[] x = Mat2Vec.mat2Vec(X, true);
		
		return x; 
	}
	
	/**
	 * transposes the matrix {@code X}
	 * <br>X_ij --> X_ji 
	 * @param X
	 * @return
	 */
	public static double[][] transpose(double[][] X){
		
		int m = X.length;
		int n = X[0].length;
		
		double[][] Y = new double[n][m];
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				Y[i][j] = X[j][i];
			}
		}
		
		return Y;		
	}
		
	/**
	 * invert matrix {@code X}
	 * @param X
	 * @return
	 */
    public static double[][] inverse(double[][] X) {
    	
    	if (X.length != X[0].length) {
    		throw new IllegalArgumentException("Only quadratic matrices can be inverted");
    	}
    	
        double[][] X_inv = new double[X.length][X.length];

        // minors and cofactors
        for (int i = 0; i < X.length; i++) {
            for (int j = 0; j < X[i].length; j++) {
                X_inv[i][j] = Math.pow(-1, i + j) * det(sub(X, i, j));
            }
        }
        // adjugate and determinant
        double det = 1.0 / det(X);
        for (int i = 0; i < X_inv.length; i++) {
            for (int j = 0; j <= i; j++) {
                double temp = X_inv[i][j];
                X_inv[i][j] = X_inv[j][i] * det;
                X_inv[j][i] = temp * det;
            }
        }

        return X_inv;
    }
	
	/**
	 * returns a submatrix by eliminating values from {@code X} in {@code row} and {@code column} 
	 * @param X
	 * @param row
	 * @param column
	 * @return
	 */
	public static double[][] sub(double[][] X, int row, int column){
		double[][] Y = new double[X.length - 1][X.length - 1];
        for (int i = 0; i < X.length; i++) {
        	for (int j = 0; i != row && j < X[i].length; j++) {
        		if (j != column) {
                    Y[i < row ? i : i - 1][j < column ? j : j - 1] = X[i][j];
        		}
        	}
        }
        return Y;
	}	
	
	/**
	 * returns the determinant of matrix {@code X}
	 * @param X
	 * @return
	 */
	public static double det(double[][] X) {
        if (X.length != X[0].length) {
            throw new IllegalArgumentException("determinant can only be computed for quadratic matrices");
        }

        if (X.length == 1 && X[0].length == 1) {
        	return X[0][0];
        }
        
        if (X.length == 2) {
            return X[0][0] * X[1][1] - X[0][1] * X[1][0];
        }
            
        double det = 0;
        for (int i = 0; i < X[0].length; i++) {
            det += Math.pow(-1, i) * X[0][i] * det(sub(X, 0, i));
        }
        return det;
    }
	
	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new columns
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[][] appendByColumns(double[][] X, double[][] Y){
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
	
	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new rows
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[][] appendByRows(double[][] X, double[][] Y){
		int n1 = X[0].length;
		int m1 = X.length;
		int n2 = Y[0].length;
		int m2 = Y.length;
		if (m1 == m2) {
			double[][] Z = new double[m1][n1 + n2];
			for (int c = 0; c < m1; c++) {
				System.arraycopy(X[c], 0, Z[c], 0, n1);
				System.arraycopy(Y[c], 0, Z[c], n1, n2);
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + ArUtils.matrixDimensionsToString(X) + " and Y" + ArUtils.matrixDimensionsToString(Y) + " do not have matching columns for appending by row!");
		}
	}
	
	private static  void checkMatrixProdDimensions(double[][] X, double[][] Y) {
		if (X[0].length != Y.length) {
            throw new IllegalArgumentException("Invalid matrix dimensions for multiplication, " + ArUtils.matrixDimensionsToString(X) + " --> " + ArUtils.matrixDimensionsToString(Y));
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
