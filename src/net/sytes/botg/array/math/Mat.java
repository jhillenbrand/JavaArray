package net.sytes.botg.array.math;

import java.util.Random;

import net.sytes.botg.array.Ar;

/**
 * Class that contains matrix to matrix operations
 * <br>if not specified in the method description, the matrices adhere the following dimension nomenclature:
 * <br>>> double[][] X = new double[r][c],
 * <br>where r is the number of rows and c the number of columns
 * @author hillenbrand
 *
 */
public class Mat {
	
	private static final double GAUSSIAN_EPSILON = 1e-20;
	
	// Suppress default constructor for noninstantiability
	private Mat() {
		throw new AssertionError(this.getClass().getSimpleName() + " cannot be instantiated");
	}
	
	/**
	 * extracts the specified {@code feature} from {@code X} and returns them in an {@code double[]} array
	 * @param X
	 * @param feature
	 * @return
	 */
	public static double[] feature(double[][] X, Feature feature) {
		if (feature != null) {
			double[] fx = new double[X.length];
			for (int i = 0; i < X.length; i++) {
				fx[i] = Vec.feature(X[i], feature);
			}
			return fx;
		} else {
			return null;
		}
	}
	
	/**
	 * extracts the specified {@code feature}s from {@code X} and returns them in an {@code double[][]} array
	 * <br>if {@code features} is specified as NULL, then all default {@code Feature}s in {@code ALL_FEATURES} are extracted
	 * <br>including features like RMS, MEAN, CREST, KURTOSIS, ... 
	 * @param X
	 * @param features
	 * @return
	 */
	public static double[][] features(double[][] X, Feature[] features) {
		if (features != null) {
			double[][] FX = new double[features.length][X.length];
			for (int f = 0; f < features.length; f++) {
				for (int i = 0; i < X.length; i++) {
					FX[f][i] = Vec.feature(X[i], features[f]);
				}
			}
			return FX;
		} else {
			return features(X, Feature.values());
		}
	}
	
	/**
	 * --------------------------------------------------------------------
	 * Matrix to Scalar
	 * --------------------------------------------------------------------
	 */

	/**
	 * return max value of all matrix elements
	 * @param X
	 * @return
	 */
	public static double max(double[][] X) {
		double xMax = Double.MIN_VALUE;
		int r = X.length;
		int c = X[0].length;
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				if (X[i][j] > xMax) {
					xMax = X[i][j];
				}
			}
		}
		return xMax;
	}
	
	/**
	 * returns the maximum values in {@code X} along dimension {@code dim}={1, 2}  
	 * @param X
	 * @param dim
	 * @return
	 */
	public static double[] max(double[][] X, int dim) {
		if (dim != 1 && dim != 2) {
			throw new IllegalArgumentException("dim must either be 1 or 2");
		}		
		Ar.checkForEqualDimensions(X);
		
		double[] max = null;
		int n = X.length;
		int m = X[0].length;
		
		if (dim == 1) {
			max = new double[n];
			for (int i = 0; i < n; i++) {
				double xMax = Double.MIN_VALUE;
				for (int j = 0; j < m; j++) {
					if (X[i][j] > xMax) {
						xMax = X[i][j];
					}
				}
				max[i] = xMax;
			}
		} else if (dim == 2){
			max = new double[X[0].length];
			for (int j = 0; j < m; j++) {
				double xMax = Double.MIN_VALUE;
				for (int i = 0; i < n; i++) {
					if (X[i][j] > xMax) {
						xMax = X[i][j];
					}
				}
				max[j] = xMax;
			}
		}
		return max;
	}
	
	/**
	 * returns the maximum values in {@code X} along dimension {@code dim}={1, 2}  
	 * @param X
	 * @param dim
	 * @return
	 */
	public static int[] max(int[][] X, int dim) {
		if (dim != 1 && dim != 2) {
			throw new IllegalArgumentException("dim must either be 1 or 2");
		}		
		Ar.checkForEqualDimensions(X);
		
		int[] max = null;
		int n = X.length;
		int m = X[0].length;
		
		if (dim == 1) {
			max = new int[n];
			for (int i = 0; i < n; i++) {
				int xMax = Integer.MIN_VALUE;
				for (int j = 0; j < m; j++) {
					if (X[i][j] > xMax) {
						xMax = X[i][j];
					}
				}
				max[i] = xMax;
			}
		} else if (dim == 2){
			max = new int[X[0].length];
			for (int j = 0; j < m; j++) {
				int xMax = Integer.MIN_VALUE;
				for (int i = 0; i < n; i++) {
					if (X[i][j] > xMax) {
						xMax = X[i][j];
					}
				}
				max[j] = xMax;
			}
		}
		return max;
	}
	
	/**
	 * returns the maximum values in {@code X} along dimension {@code dim}={1, 2}  
	 * @param X
	 * @param dim
	 * @return
	 */
	public static float[] max(float[][] X, int dim) {
		if (dim != 1 && dim != 2) {
			throw new IllegalArgumentException("dim must either be 1 or 2");
		}		
		Ar.checkForEqualDimensions(X);
		
		float[] max = null;
		int n = X.length;
		int m = X[0].length;
		
		if (dim == 1) {
			max = new float[n];
			for (int i = 0; i < n; i++) {
				float xMax = Float.MIN_VALUE;
				for (int j = 0; j < m; j++) {
					if (X[i][j] > xMax) {
						xMax = X[i][j];
					}
				}
				max[i] = xMax;
			}
		} else if (dim == 2){
			max = new float[X[0].length];
			for (int j = 0; j < m; j++) {
				float xMax = Float.MIN_VALUE;
				for (int i = 0; i < n; i++) {
					if (X[i][j] > xMax) {
						xMax = X[i][j];
					}
				}
				max[j] = xMax;
			}
		}
		return max;
	}	
	
	/**
	 * rms values of {@code X} along dimension {@code dim}
	 * @param X
	 * @param dim
	 * @return
	 */
	public static double[] rms(double[][] X, int dim) {
		if (dim != 1 && dim != 2) {
			throw new IllegalArgumentException("dim must either be 1 or 2");
		}		
		Ar.checkForEqualDimensions(X);
		
		double[] rms = null;
		int n = X.length;
		int m = X[0].length;
		
		if (dim == 1) {
			rms = new double[n];
			for (int i = 0; i < n; i++) {
				double rmsSum = 0;
				for (int j = 0; j < m; j++) {
					rmsSum = rmsSum + X[i][j] * X[i][j];
				}
				rms[i] = Math.sqrt(rmsSum / m);
			}
		} else if (dim == 2) {
			rms = new double[m];
			for (int j = 0; j < m; j++) {
				double rmsSum = 0;
				for (int i = 0; i < n; i++) {
					rmsSum = rmsSum + X[i][j] * X[i][j];
				}
				rms[j] = Math.sqrt(rmsSum / n);
			}
		}
		return rms;
	}
	
	/**
	 * rms values of {@code X} along dimension {@code dim}
	 * @param X
	 * @param dim
	 * @return
	 */
	public static float[] rms(float[][] X, int dim) {
		if (dim != 1 && dim != 2) {
			throw new IllegalArgumentException("dim must either be 1 or 2");
		}		
		Ar.checkForEqualDimensions(X);
		
		float[] rms = null;
		int n = X.length;
		int m = X[0].length;
		
		if (dim == 1) {
			rms = new float[n];
			for (int i = 0; i < n; i++) {
				float rmsSum = 0;
				for (int j = 0; j < m; j++) {
					rmsSum = rmsSum + X[i][j] * X[i][j];
				}
				rms[i] = (float) Math.sqrt(rmsSum / m);
			}
		} else if (dim == 2) {
			rms = new float[m];
			for (int j = 0; j < m; j++) {
				float rmsSum = 0;
				for (int i = 0; i < n; i++) {
					rmsSum = rmsSum + X[i][j] * X[i][j];
				}
				rms[j] = (float) Math.sqrt(rmsSum / n);
			}
		}
		return rms;
	}
	
	/**
	 * rms values of {@code X} along dimension {@code dim}
	 * @param X
	 * @param dim
	 * @return
	 */
	public static int[] rms(int[][] X, int dim) {
		if (dim != 1 && dim != 2) {
			throw new IllegalArgumentException("dim must either be 1 or 2");
		}		
		Ar.checkForEqualDimensions(X);
		
		int[] rms = null;
		int n = X.length;
		int m = X[0].length;
		
		if (dim == 1) {
			rms = new int[n];
			for (int i = 0; i < n; i++) {
				int rmsSum = 0;
				for (int j = 0; j < m; j++) {
					rmsSum = rmsSum + X[i][j] * X[i][j];
				}
				rms[i] = (int) Math.sqrt(rmsSum / m);
			}
		} else if (dim == 2) {
			rms = new int[m];
			for (int j = 0; j < m; j++) {
				int rmsSum = 0;
				for (int i = 0; i < n; i++) {
					rmsSum = rmsSum + X[i][j] * X[i][j];
				}
				rms[j] = (int) Math.sqrt(rmsSum / n);
			}
		}
		return rms;
	}

	/**
	 * returns the minimum values in {@code X} along dimension {@code dim}={1, 2}  
	 * @param X
	 * @param dim
	 * @return
	 */
	public static double[] min(double[][] X, int dim) {
		if (dim != 1 && dim != 2) {
			throw new IllegalArgumentException("dim must either be 1 or 2");
		}		
		Ar.checkForEqualDimensions(X);
		
		double[] min = null;
		int n = X.length;
		int m = X[0].length;
		
		if (dim == 1) {
			min = new double[n];
			for (int i = 0; i < n; i++) {
				double xMin = Double.MAX_VALUE;
				for (int j = 0; j < m; j++) {
					if (X[i][j] < xMin) {
						xMin = X[i][j];
					}
				}
				min[i] = xMin;
			}
		} else if (dim == 2){
			min = new double[X[0].length];
			for (int j = 0; j < m; j++) {
				double xMin = Double.MAX_VALUE;
				for (int i = 0; i < n; i++) {
					if (X[i][j] < xMin) {
						xMin = X[i][j];
					}
				}
				min[j] = xMin;
			}
		}
		return min;
	}

	/**
	 * returns the minimum values in {@code X} along dimension {@code dim}={1, 2}  
	 * @param X
	 * @param dim
	 * @return
	 */
	public static float[] min(float[][] X, int dim) {
		if (dim != 1 && dim != 2) {
			throw new IllegalArgumentException("dim must either be 1 or 2");
		}		
		Ar.checkForEqualDimensions(X);
		
		float[] min = null;
		int n = X.length;
		int m = X[0].length;
		
		if (dim == 1) {
			min = new float[n];
			for (int i = 0; i < n; i++) {
				float xMin = Float.MAX_VALUE;
				for (int j = 0; j < m; j++) {
					if (X[i][j] < xMin) {
						xMin = X[i][j];
					}
				}
				min[i] = xMin;
			}
		} else if (dim == 2){
			min = new float[X[0].length];
			for (int j = 0; j < m; j++) {
				float xMin = Float.MAX_VALUE;
				for (int i = 0; i < n; i++) {
					if (X[i][j] < xMin) {
						xMin = X[i][j];
					}
				}
				min[j] = xMin;
			}
		}
		return min;
	}

	/**
	 * returns the minimum values in {@code X} along dimension {@code dim}={1, 2}  
	 * @param X
	 * @param dim
	 * @return
	 */
	public static int[] min(int[][] X, int dim) {
		if (dim != 1 && dim != 2) {
			throw new IllegalArgumentException("dim must either be 1 or 2");
		}		
		Ar.checkForEqualDimensions(X);
		
		int[] min = null;
		int n = X.length;
		int m = X[0].length;
		
		if (dim == 1) {
			min = new int[n];
			for (int i = 0; i < n; i++) {
				int xMin = Integer.MAX_VALUE;
				for (int j = 0; j < m; j++) {
					if (X[i][j] < xMin) {
						xMin = X[i][j];
					}
				}
				min[i] = xMin;
			}
		} else if (dim == 2){
			min = new int[X[0].length];
			for (int j = 0; j < m; j++) {
				int xMin = Integer.MAX_VALUE;
				for (int i = 0; i < n; i++) {
					if (X[i][j] < xMin) {
						xMin = X[i][j];
					}
				}
				min[j] = xMin;
			}
		}
		return min;
	}
	
	/**
	 * return min value of all matrix elements
	 * @param X
	 * @return
	 */
	public static double min(double[][] X) {
		Ar.checkForNull(X);
		double xMin = Double.MAX_VALUE;
		int r = X.length;
		int c = X[0].length;
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				if (X[i][j] < xMin) {
					xMin = X[i][j];
				}
			}
		}
		return xMin;
	}
	
	public static double mean(double[][] X) {
		Ar.checkForNull(X);
		double sum = 0.0;
		int r = X.length;
		int c = X[0].length;
		int s = r * c;
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				sum = sum + X[i][j];
			}
		}
		return sum / s;
	}
	
	public static double variance(double[][] X) {
		Ar.checkForNull(X);
		int r = X.length;
		int c = X[0].length;
		double mean = mean(X);
		double variance = 0.0;
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				variance = variance + Math.pow(X[i][j] - mean, 2);
			}
		}
		return variance / (r * c);
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
            det += Math.pow(-1, i) * X[0][i] * det(subByElimination(X, 0, i));
        }
        return det;
    } 
	
	/**
	 * return min and max as array elements of all matrix elements,
	 * <br>Example:
	 * <br>
	 * <br>double[] mm = minMax(x);
	 * <br>min = mm[0];
	 * <br>max = mm[1]; 
	 * @param X
	 * @return
	 */
	public static double[] minMax(double[][] X) {
		Ar.checkForNull(X);
		double xMin = Double.MAX_VALUE;
		double xMax = Double.MIN_VALUE;
		int r = X.length;
		int c = X[0].length;
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				if (X[i][j] < xMin) {
					xMin = X[i][j];
				}
				if (X[i][j] > xMax) {
					xMax = X[i][j];
				}
			}
		}
		return new double[]{xMin, xMax};
	}
	
	/**
	 * --------------------------------------------------------------------
	 * Matrix to Vector
	 * --------------------------------------------------------------------
	 */

	/**
	 * solve the linear equation system y<sub>j</sub> = A<sub>ij</sub> x<sub>i</sub> for x<sub>i</sub>
	 * <br>using Gaussian elimination with partial pivoting
	 * <br>adapted from <a href="https://introcs.cs.princeton.edu/java/95linear/GaussianElimination.java.html">LINK</a>
	 * @param A
	 * @param y
	 * @return
	 */
	public static double[] linSolveGaussian(double[][] A, double[] y) {
		if (A.length != y.length) {
			throw new IllegalArgumentException("First matrix dimension of A (" + A.length + ") and length of vector y (" + y.length + ") must be equal");
		}
		
		int n = y.length;
		
		for (int p = 0; p < n; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            double[] temp = A[p];
            A[p] = A[max];
            A[max] = temp;
            double t = y[p];
            y[p] = y[max];
            y[max] = t;

            // singular or nearly singular
            if (Math.abs(A[p][p]) <= GAUSSIAN_EPSILON) {
                throw new IllegalArgumentException("Matrix is singular or nearly singular");
            }

            // pivot within A and b
            for (int i = p + 1; i < n; i++) {
                double alpha = A[i][p] / A[p][p];
                y[i] -= alpha * y[p];
                for (int j = p; j < n; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }
		
		// back substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (y[i] - sum) / A[i][i];
        }
        return x;
	}
	
	/**
	 * solve the linear equation system y<sub>j</sub> = A<sub>ij</sub> x<sub>i</sub> for x<sub>i</sub>
	 * <br>NOTE: only works in moderate time for small problems j < 10
	 * @param A
	 * @param y
	 * @return
	 */
	public static double[] linSolveInverse(double[][] A, double[] y) {
		
		if (A.length != y.length) {
			throw new IllegalArgumentException("First matrix dimension of A (" + A.length + ") and length of vector y (" + y.length + ") must be equal");
		}
		
		int n = A.length;
		
		// put y-vector into matrix Y for matrix multiplication
		double[][] Y = Vec.matrix(y, n, 1, false);
		
		double[][] X = product(Mat.inverse(A), Y); 
		
		// reshape the result matrix back to vector
		double[] x = Mat.vector(X, true);
		
		return x; 
	}	
	
	/**
	 * TODO NOT WORKING
	 * @param A
	 * @param b
	 * @return
	 */
	public static double[] linSolveLU(double[][] A, double[] b) {
		Ar.checkForNull(A);
		Ar.checkForNull(b);
		checkMatrixProdDimensions(A, b);
		int n = A.length;
        // decomposition of matrix
        double[][] LU = new double[n][n];
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += LU[i][k] * LU[k][j];
                }
                LU[i][j] = A[i][j] - sum;
            }
            for (int j = i + 1; j < n; j++){
                sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += LU[j][k] * LU[k][i];
                }
                LU[j][i] = (1 / LU[i][i]) * (A[j][i] - sum);
            }
        }

        // lu = L+U-I
        // find solution of Ly = b
         double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            sum = 0;
            for (int k = 0; k < i; k++) {
                sum += LU[i][k] * y[k];
            }
            y[i] = b[i] - sum;
        }
        // find solution of Ux = y
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--){
            sum = 0;
            for (int k = i + 1; k < n; k++) {
                sum += LU[i][k] * x[k];
            }
            x[i] = (1 / LU[i][i]) * (y[i] - sum);
        }
        return x;
    }
	
	/**
	 * matrix-vector multiplication
	 * @param A
	 * @param x
	 * @return
	 */
	public static double[] product(double[][] A, double[] x) {
		double[][] X = Vec.matrix(x, 1);
		double[][] X_T = Mat.transpose(X);
		double[][] Y = product(A, X_T);
		double[][] Y_T = transpose(Y);
		return Y_T[0];
	}
	
	/**
	 * --------------------------------------------------------------------
	 * Matrix to Matrix
	 * --------------------------------------------------------------------
	 */

	/**
	 * transposes the matrix {@code X}
	 * <br>X_ij --> X_ji 
	 * @param X
	 * @return
	 */
	public static double[][] transpose(double[][] X){
		
		int r = X.length;
		int c = X[0].length;
		
		double[][] Y = new double[c][r];
		
		for (int i = 0; i < c; i++) {
			for (int j = 0; j < r; j++) {
				Y[i][j] = X[j][i];
			}
		}
		
		return Y;		
	}
	
	/**
	 * transposes the matrix {@code X}
	 * <br>X_ij --> X_ji 
	 * @param X
	 * @return
	 */
	public static float[][] transpose(float[][] X){
		
		int r = X.length;
		int c = X[0].length;
		
		float[][] Y = new float[c][r];
		
		for (int i = 0; i < c; i++) {
			for (int j = 0; j < r; j++) {
				Y[i][j] = X[j][i];
			}
		}
		
		return Y;		
	}
	
	/**
	 * transposes the matrix {@code X}
	 * <br>X_ij --> X_ji 
	 * @param X
	 * @return
	 */
	public static int[][] transpose(int[][] X){
		
		int r = X.length;
		int c = X[0].length;
		
		int[][] Y = new int[c][r];
		
		for (int i = 0; i < c; i++) {
			for (int j = 0; j < r; j++) {
				Y[i][j] = X[j][i];
			}
		}
		
		return Y;		
	}
	
	/**
	 * transposes the matrix {@code X}
	 * <br>X_ij --> X_ji 
	 * @param X
	 * @return
	 */
	public static long[][] transpose(long[][] X){
		
		int r = X.length;
		int c = X[0].length;
		
		long[][] Y = new long[c][r];
		
		for (int i = 0; i < c; i++) {
			for (int j = 0; j < r; j++) {
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
                X_inv[i][j] = Math.pow(-1, i + j) * det(subByElimination(X, i, j));
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
     * TODO NOT WORKING
     * returns an approximate form of the inverse to matrix {@code X} using the Neumann-Series
     * @param X matrix to be inverted
     * @param k number of series elements
     * @return
     */
    public static double[][] inverseNumerical(double[][] X, int k){
    	if (!isSquare(X)) {
    		throw new IllegalArgumentException("matrix X must be square");
    	}
    	int n = X.length;
    	double[][] X_inv = new double[n][n];
    	double[][] I = unitMatrix(n);
    	double[][] T = minus(I, X);
    	for (int i = 0; i < k; i++) {
    		X_inv = plus(X_inv, power(T, i));
    	}
    	return X_inv;
    }
	
	/**
	 * returns a submatrix by eliminating elements from {@code X} in {@code row} and {@code column} 
	 * @param X
	 * @param row
	 * @param column
	 * @return
	 */
	public static double[][] subByElimination(double[][] X, int row, int column){
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
		Ar.checkForNull(points);
		Ar.checkForEmpty(points);
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
				double d = Vec.distance(points[i], points[j]);
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
		Ar.checkForNull(XYZ);
		Ar.checkForEmpty(XYZ);
		int r = XYZ.length;
		int c = XYZ[0].length;
		double[][] XYZ_n = new double[r][c];
		for (int j = 0; j < c; j++) {
			double squaredSum = 0.0;
			for (int i = 0; i < r; i++) {
				squaredSum = squaredSum + XYZ[i][j] * XYZ[i][j]; 
			}
			double vecLen = Math.sqrt(squaredSum);
			for (int i = 0; i < r; i++) {
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
		int r = X.length;
		int c = X[0].length;		
		double xMin = min(X);
		double xMax = max(X);
		double[][] xn = new double[r][c];
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				xn[i][j] = (X[i][j] - xMin) / (xMax - xMin);
			}
		}
		return xn;
	} 
	
	/**
	 * normalization of matrix {@code X} with zscore
	 * @param X
	 * @return
	 */
	public static double[][] zscore(double[][] X) {
		Ar.checkForEmpty(X);
		double mean = mean(X);
		double sigma = variance(X);
		int r = X.length;
		int c = X[0].length;
		double[][] Z = new double[r][c];
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				Z[i][j] = (X[i][j] - mean) / sigma;
			}
		}
		return Z;
	}
	
	/**
	 * matrix addition  {@code X} + {@code Y}
	 * <br>Z<sub>ij</sub> = X<sub>ij</sub> + Y<sub>ij</sub>
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[][] plus(double[][] X, double[][] Y){
		Ar.checkForNull(X, Y);
		Ar.checkForEmpty(X, Y);
		Ar.checkForMatchingDimensions(X, Y);
		int r = X.length;
		int c = X[0].length;
		double[][] Z = new double[r][c];
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				Z[i][j] = X[i][j] + Y[i][j];
			}
		}
		return Z;
	}
	
	/**
	 * matrix subtraction {@code X} - {@code Y}
	 * <br>Z<sub>ij</sub> = X<sub>ij</sub> - Y<sub>ij</sub>
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[][] minus(double[][] X, double[][] Y){
		Ar.checkForNull(X, Y);
		Ar.checkForEmpty(X, Y);
		Ar.checkForMatchingDimensions(X, Y);
		int r = X.length;
		int c = X[0].length;
		double[][] Z = new double[r][c];
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				Z[i][j] = X[i][j] - Y[i][j];
			}
		}
		return Z;
	}
	
	/**
	 * vector subraction x<sub>i</sub> from each matrix element A<sub>ij</sub>
	 * @param A
	 * @param x
	 * @return
	 */
	public static double[][] minus(double[][] A, double[] x){
		if (A.length != x.length) {
			throw new IllegalArgumentException("First dimension of A[" + A.length + "][" + A[0].length + "] does not match dimension of x[" + x.length + "].");
		}
		int r = A.length;
		int c = A[0].length;
		double[][] B = new double[r][c];
		for (int i = 0; i < r; i++) {
			double d = x[i];
			for (int j = 0; j < c; j++) {
				B[i][j] = A[i][j] - d;
			}
		}
		return B;
	}
	
	/**
	 * returns the matrix multiplication of {@code X} and {@code Y}
	 * <br>Z<sub>ij</sub> = X<sub>ik</sub>Y<sub>kj</sub>
	 * <br>using naive approach of 3 nested loops i -> j -> k
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[][] product3(double[][] X, double[][] Y){
		Ar.checkForNull(X, Y);
		checkMatrixProdDimensions(X, Y);
		int rX = X.length;
		int cX = X[0].length;
		int cY = Y[0].length;
        double[][] Z = new double[rX][cY];
        for (int i = 0; i < rX; i++) {
            for (int j = 0; j < cY; j++) {
                double sum = 0;
                for (int k = 0; k < cX; k++) {
                    sum += X[i][k] * Y[k][j];
                }
                Z[i][j] = sum;
            }
        }
        return Z;
	}
	
	/**
	 * returns the matrix multiplication of {@code X} and {@code Y}
	 * <br>Z<sub>ij</sub> = X<sub>ik</sub>Y<sub>kj</sub>
	 * <br>using JAMA's algorithm adapted from <a href="http://www.ii.uib.no/~geirg/NIK2002.pdf">LINK</a>
	 * @param X 1st matrix
	 * @param Y 2nd matrix
	 * @return Z<sub>ij</sub>, 2D array
	 */
	public static double[][] product2(double[][] X, double[][] Y){
		Ar.checkForNull(X, Y);
		checkMatrixProdDimensions(X, Y);
		int rX = X.length;
		int cX = X[0].length;
		int rY = Y.length;
		int cY = Y[0].length;
		double[][] Z = new double[rX][cY];
		double[] Ycolj = new double[cX]; // for caching 
		for (int j = 0; j < cY; j++) {
			for (int k = 0; k < cX; k++) {
				Ycolj[k] = Y[k][j];
			}
			for (int i = 0; i < rX; i++) {
				double[] Xrowi = X[i];
				double sum = 0;
				for (int k = 0; k < cX; k++) {
					sum += Xrowi[k] * Ycolj[k];
				}
				Z[i][j] = sum;
			}
		}		
		return Z;
	}
	
	/**
	 * returns the matrix multiplication of {@code X} and {@code Y}
	 * <br>Z<sub>ij</sub> = X<sub>ik</sub>Y<sub>kj</sub>
	 * <br>using pure row oriented algorithm (i,k,j) adapted from <a href="http://www.ii.uib.no/~geirg/NIK2002.pdf">LINK</a>
	 * @param X 1st matrix
	 * @param Y 2nd matrix
	 * @return Z<sub>ij</sub>, 2D array
	 */
	public static double[][] product(double[][] X, double[][] Y) {
		Ar.checkForNull(X, Y);
		checkMatrixProdDimensions(X, Y);
		int rX = X.length;
		int cX = X[0].length;
		int rY = Y.length;
		int cY = Y[0].length;
		double[][] Z = new double[rX][cY];
		double[] Xrowi, Yrowi, Zrowi;
		int i = 0, j = 0, k = 0;
		double a = 0.0;
		for (i = 0; i < rX; i++) {
			Xrowi = X[i];
			Zrowi = Z[i];
			for (k = 0; k < rY; k++){
				Yrowi = Y[k];
				a = Xrowi[k];
				for (j = cY; --j >= 0;) {
					Zrowi[j] += a * Yrowi[j];
				}
			}
		}
		return Z;
	}
	
	public static double[][] tiledProduct(double[][] X, double[][] Y, int tileSize) {
		Ar.checkForNull(X, Y);
		checkMatrixProdDimensions(X, Y);
		int rX = X.length;
		int cX = X[0].length;
		int cY = Y[0].length;
        double[][] Z = new double[rX][cY];
        for (int ih = 0; ih < rX; ih += tileSize) {
            for (int jh = 0; jh < cY; jh += tileSize) {
                for (int kh = 0; kh < cX; kh += tileSize) {
                	for (int il = 0; il < tileSize; ++il) {
                		for (int kl = 0; kl < tileSize; ++kl) {
                			for (int jl = 0; jl < tileSize; ++jl) {
                				Z[ih+il][jh+jl] += X[ih+il][kh+kl] * Y[kh+kl][jh+jl];
                			}
                		}
                	}
                }
            }
        }
        return Z;
	}
	
	
	/**
	 * returns the matrix multiplication of {@code X} and {@code Y}
	 * <br>Z<sub>ij</sub> = X<sub>ik</sub>Y<sub>kj</sub>
	 * <br>using naive approach of 3 nested loops i -> k -> j
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[][] product4(double[][] X, double[][] Y){
		Ar.checkForNull(X, Y);
		checkMatrixProdDimensions(X, Y);
		int rX = X.length;
		int cX = X[0].length;
		int cY = Y[0].length;
        double[][] Z = new double[rX][cY];
        for (int i = 0; i < rX; i++) {
            for (int k = 0; k < cX; k++) {
                for (int j = 0; j < cY; j++) {
                	Z[i][j] += X[i][k] * Y[k][j];
                }
            }
        }
        return Z;
	}
	
	
	/**
	 * multiplies a scalar {@code c} with every element in {@code X}
	 * @param X
	 * @param x
	 * @return
	 */
	public static double[][] product(double[][] X, double c){
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		Ar.checkForEqualDimensions(X);
		double[][] Y = new double[X.length][X[0].length];
		for (int i = 0; i < X.length; i++) {
			for (int j = 0; j < X[0].length; j++) {
				Y[i][j] = X[i][j] * c;
			}
		}
		return Y;
	}
	
	/**
	 * multiplies a scalar {@code c} with every element in {@code X}
	 * this changes the matrix that was passed;
	 * @param X
	 * @param x
	 * @return
	 */
	public static void product2(double[][] X, double c){
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		Ar.checkForEqualDimensions(X);
		for (int i = 0; i < X.length; i++) {
			for (int j = 0; j < X[0].length; j++) {
				X[i][j] = X[i][j] * c;
			}
		}
	}	
	
	
	public static double[][] square(double[][] X){
		if (isSquare(X)) {
			return product(X, X);
		} else {
			throw new IllegalArgumentException("matrix must be square");
		}
	}
	
	/**
	 * returns the exponentiation of matrix {@code X} by {@code k} using exponentiation by squaring
	 * @param X
	 * @param k
	 * @return
	 */
	public static double[][] power(double[][] X, int k){
		if (isSquare(X)) {
			if (k < 0) {
				throw new IllegalArgumentException("non-positive exponents are not defined");
			} else if (k == 0 ) {
				return unitMatrix(X.length);
			} else if(k == 1) {
				return Ar.copy(X);
			} else {
				if (k % 2 == 0) {
					double[][] T = product(X, X);
					return power(T, k / 2);
				} else {
					double[][] T = product(X, product(X, X));
					return power(T, (k - 1) / 2);
				}
			}
		} else {
			throw new IllegalArgumentException("matrix must be square");
		}
	}
		
	/**
	 * --------------------------------------------------------------------
	 * Matrix Generation Methods
	 * -------------------------------------------------------------------- 
	 */
	
	/**
	 * return an array with zeros of size [r, c]
	 * @param r number of rows
	 * @param c number of columns
	 * @return
	 */
	public static double[][] zeros(int r, int c) {
		return new double[r][c];
	}
	
	/**
	 * return an float array with zeros of size [r, c]
	 * @param r number of rows
	 * @param c number of columns
	 * @return
	 */
	public static float[][] zerosF(int r, int c) {
		return new float[r][c];
	}
		
	/**
	 * returns an 2D array with 1's, where {@code r} is the number of rows and {@code c} the number of columns
	 * @param r number of rows
	 * @param c number of columns
	 * @return
	 */
	public static double[][] ones(int r, int c){
		double[][] ar = new double[r][c];
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				ar[i][j] = 1.0;
			}
		}
		return ar;
	}
	
	/**
	 * returns an 2D float array with 1's, where {@code r} is the number of rows and {@code c} the number of columns
	 * @param r
	 * @param c
	 * @return
	 */
	public static float[][] onesF(int r, int c){
		float[][] ar = new float[r][c];
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				ar[i][j] = 1.0f;
			}
		}
		return ar;
	}
	
	/**
	 * returns an 2D int array with 1's, where {@code r} is the number of rows and {@code c} the number of columns
	 * @param r
	 * @param c
	 * @return
	 */
	public static int[][] onesI(int r, int c){
		int[][] ar = new int[r][c];
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				ar[i][j] = 1;
			}
		}
		return ar;
	}
	
	
	/**
	 * returns a squared unit matrix of size {@code n}
	 * <br> Example for n = 3
	 * <br>| 1, 0, 0 |
	 * <br>| 0, 1, 0 | 
	 * <br>| 0, 0, 1 |
	 * @param n
	 * @return
	 */
	public static double[][] unitMatrix(int n){
		double[][] ar = zeros(n, n);
		for (int i = 0; i < n; i++) {
			ar[i][i] = 1;
		}
		return ar;
	}
	
	/**
	 * returns a unit square float matrix of size {@code n}
	 * <br> Example for n = 3
	 * <br>| 1f, 0,  0 |
	 * <br>| 0, 1f,  0 | 
	 * <br>| 0,  0, 1f |
	 * @param n
	 * @return
	 */
	public static float[][] unitMatrixF(int n){
		float[][] ar = zerosF(n, n);
		for (int i = 0; i < n; i++) {
			ar[i][i] = 1;
		}
		return ar;
	}	
	
	/**
	 * returns a 2D matrix with random values in [0, 1] of size [{@code r},{@code c}]
	 * @param r number of rows
	 * @param c number of columns
	 * @return
	 */
	public static double[][] rand(int r, int c){
		Random rr = new Random();
		double[][] ar = new double[r][c];
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				ar[i][j] = rr.nextDouble();
			}
		}
		return ar;
	}
		
	/**
	 * returns a 2D matrix with random values from all 2^32 possible int values of size [{@code r},{@code c}]
	 * @param r number of rows
	 * @param c number of columns
	 * @return
	 */
	public static int[][] randI(int r, int c){
		Random rr = new Random();
		int[][] ar = new int[r][c];
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				ar[i][j] = rr.nextInt();
			}
		}
		return ar;
	}
	
	
	/**
	 * returns a matrix X � R<sup>m x n</sup> with elements increasing from 0 to k = {@code r} * {@code c}
	 * @param r number of rows
	 * @param c number of columns
	 * @return
	 */
	public static double[][] incrementMat(int r, int c) {
		int k = r * c;
		double[][] X = new double[r][c];		
		for (int i = 0; i < k; i++) {
			int rem = i / c;
			int n = i % c;
			X[rem][n] = i;
		}
		return X;
	}
		
	/**
	 * creates a meshgrid of the two input vectors
	 * @param x
	 * @param y
	 * @return
	 */
	public static double[][][] meshgrid(double[] x, double[] y){
		
		int n = x.length;
		int m = y.length;
		
		double[][] X = new double[n][m];
		double[][] Y = new double[n][m];
		
		for (int i = 0; i < n; i++) {			
			for (int j = 0; j < m; j++) {
				X[i][j] = x[i];
				Y[i][j] = y[j];
			}			
		}
		
		double[][][] XY = new double[2][][];
		XY[0] = X;
		XY[1] = Y;
		
		return XY;
		
	}
	
	/**
	 * returns 2D matrix with NaN values of size [{@code r} x {@code c}]
	 * @param r
	 * @param c
	 * @return
	 */
	public static double[][] nan(int r, int c) {
		double[][] mat = new double[r][c];
		for (int i = 0; i < r; i++) {
			mat[i] = Vec.nan(c);
		}
		return mat;
	}
	
	/**
	 * TODO
	 * creates the Levi-Civita Matrix &epsilon;<sub>ijk...</sub> based on dimension {@code d}
	 * @param d
	 * @return
	 */
	public static double[][] leviCivita(int d){
		if (d == 1 || d > 3) {
			throw new IllegalArgumentException("");
		}
		return null;
	}
	
	
	
	/**
	 * --------------------------------------------------------------------
	 * Matrix Utility Methods
	 * --------------------------------------------------------------------
	 */

	
	/**
	 * returns the number of elements in {@code X}
	 * @param X
	 * @return
	 */
	public static int numel(double[][] X) {
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		int n = 0;
		for (int i = 0; i < X.length; i++) {
			n = n + X[i].length;
		}
		return n;
	}
	
	/**
     * checks whether the matrix {@code X} is square
     * <br>{@code X} &#8712; R<sup>n x n</sup> 
     * @param X
     * @return
     */
	public static boolean isSquare(double[][] X) {
		if (X.length != X[0].length) {
			return false; 
		} else {
			return true;
		}
	}
	
	/**
	 * returns the matrix {@code X} as vector, where all matrix elements are listed column by column
	 * @param X
	 * @return
	 */
	public static double[] vector(double[][] X) {
		return vector(X, true);
	}
	
	/**
	 * 
	 * @param X
	 * @return
	 */
	public static double[] vector(double[][] X, boolean columnByColumn) {
		Ar.checkForNull(X);
		Ar.checkForEmpty(X);
		int r = X.length;
		int c = X[0].length;
		double[] x = new double[r * c];
		if (columnByColumn) {
			int n = 0;
			for (int j = 0; j < c; j++) {
				for (int i = 0; i < r; i++) {
					x[n] = X[i][j];
					++n;
				}
			}
		} else {
			int n = 0;
			for (int i = 0; i < r; i++) {
				for (int j = 0; j < c; j++) {
					x[n] = X[i][j];
					++n;
				}
			}
		}
		return x;
	}
	
	/**
	 * returns the trajectory matrix (Hankel matrix) X for {@code x} and window length {@code L} 
	 * @see <a href="https://en.wikipedia.org/wiki/Singular_spectrum_analysis">wikipedia</a>
	 * @param x
	 * @param L
	 * @return
	 */
	public static double[][] trajectoryMatrix(double[] x, int L){
		int N = x.length;
		int K = N - (L - 1);
		double[][] X = new double[L][K];
		int s = 0;
		for (int k = 0; k < K; k++) {			
			for (int l = 0; l < L; l++) {
				X[l][k] = x[l + s];
			}
			++s;
		}
		return X;
	}
	
	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new rows
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[][] appendRows(double[][] X, double[][] Y){
		int c1 = X[0].length;
		int r1 = X.length;
		int c2 = Y[0].length;
		int r2 = Y.length;
		if (c1 == c2) {
			double[][] Z = new double[r1 + r2][c1];
			int n;
			for (n = 0; n < r1; n++) {
				Z[n] = X[n];
			}
			int n2 = n;
			for (n = 0; n < r2; n++) {
				Z[n2] = Y[n];
				n2 = n2 + n;
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + matrixDimensionsToString(X) + " and Y" + matrixDimensionsToString(Y) + " do not have matching columns for appending by rows!");
		}
	}

	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new rows
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static char[][] appendRows(char[][] X, char[][] Y){
		int c1 = X[0].length;
		int r1 = X.length;
		int c2 = Y[0].length;
		int r2 = Y.length;
		if (c1 == c2) {
			char[][] Z = new char[r1 + r2][c1];
			int n;
			for (n = 0; n < r1; n++) {
				Z[n] = X[n];
			}
			int n2 = n;
			for (n = 0; n < r2; n++) {
				Z[n2] = Y[n];
				n2 = n2 + n;
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + matrixDimensionsToString(X) + " and Y" + matrixDimensionsToString(Y) + " do not have matching columsn for appending by rows!");
		}
	}
	
	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new rows
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static float[][] appendRows(float[][] X, float[][] Y){
		int c1 = X[0].length;
		int r1 = X.length;
		int c2 = Y[0].length;
		int r2 = Y.length;
		if (c1 == c2) {
			float[][] Z = new float[r1 + r2][c1];
			int n;
			for (n = 0; n < r1; n++) {
				Z[n] = X[n];
			}
			int n2 = n;
			for (n = 0; n < r2; n++) {
				Z[n2] = Y[n];
				n2 = n2 + n;
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + matrixDimensionsToString(X) + " and Y" + matrixDimensionsToString(Y) + " do not have matching columsn for appending by rows!");
		}
	}
	
	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new rows
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static long[][] appendRows(long[][] X, long[][] Y){
		int c1 = X[0].length;
		int r1 = X.length;
		int c2 = Y[0].length;
		int r2 = Y.length;
		if (c1 == c2) {
			long[][] Z = new long[r1 + r2][c1];
			int n;
			for (n = 0; n < r1; n++) {
				Z[n] = X[n];
			}
			int n2 = n;
			for (n = 0; n < r2; n++) {
				Z[n2] = Y[n];
				n2 = n2 + n;
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + matrixDimensionsToString(X) + " and Y" + matrixDimensionsToString(Y) + " do not have matching columsn for appending by rows!");
		}
	}
	
	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new rows
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static int[][] appendRows(int[][] X, int[][] Y){
		int c1 = X[0].length;
		int r1 = X.length;
		int c2 = Y[0].length;
		int r2 = Y.length;
		if (c1 == c2) {
			int[][] Z = new int[r1 + r2][c1];
			int n;
			for (n = 0; n < r1; n++) {
				Z[n] = X[n];
			}
			int n2 = n;
			for (n = 0; n < r2; n++) {
				Z[n2] = Y[n];
				n2 = n2 + n;
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + matrixDimensionsToString(X) + " and Y" + matrixDimensionsToString(Y) + " do not have matching columsn for appending by rows!");
		}
	}	

	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new rows
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static String[][] appendRows(String[][] X, String[][] Y){
		int c1 = X[0].length;
		int r1 = X.length;
		int c2 = Y[0].length;
		int r2 = Y.length;
		if (c1 == c2) {
			String[][] Z = new String[r1 + r2][c1];
			int n;
			for (n = 0; n < r1; n++) {
				Z[n] = X[n];
			}
			int n2 = n;
			for (n = 0; n < r2; n++) {
				Z[n2] = Y[n];
				n2 = n2 + n;
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + matrixDimensionsToString(X) + " and Y" + matrixDimensionsToString(Y) + " do not have matching columsn for appending by rows!");
		}
	}

	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new rows
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static byte[][] appendRows(byte[][] X, byte[][] Y){
		int c1 = X[0].length;
		int r1 = X.length;
		int c2 = Y[0].length;
		int r2 = Y.length;
		if (c1 == c2) {
			byte[][] Z = new byte[r1 + r2][c1];
			int n;
			for (n = 0; n < r1; n++) {
				Z[n] = X[n];
			}
			int n2 = n;
			for (n = 0; n < r2; n++) {
				Z[n2] = Y[n];
				n2 = n2 + n;
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + matrixDimensionsToString(X) + " and Y" + matrixDimensionsToString(Y) + " do not have matching columsn for appending by rows!");
		}
	}		

	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new rows
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static boolean[][] appendRows(boolean[][] X, boolean[][] Y){
		int c1 = X[0].length;
		int r1 = X.length;
		int c2 = Y[0].length;
		int r2 = Y.length;
		if (c1 == c2) {
			boolean[][] Z = new boolean[r1 + r2][c1];
			int n;
			for (n = 0; n < r1; n++) {
				Z[n] = X[n];
			}
			int n2 = n;
			for (n = 0; n < r2; n++) {
				Z[n2] = Y[n];
				n2 = n2 + n;
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + matrixDimensionsToString(X) + " and Y" + matrixDimensionsToString(Y) + " do not have matching columsn for appending by rows!");
		}
	}	
	
	/**
	 * creates a new matrix by appending the vector {@code x} to matrix {@code X} as new row
	 * @param X
	 * @param x
	 * @return
	 */
	public static double[][] appendRows(double[][] X, double[] x){
		double[][] Y = new double[X.length + 1][];
		for (int i = 0; i < X.length; i++) {
			Y[i] = X[i];
		}
		Y[Y.length - 1] = x; 
		return Y;
	}
	
	/**
	 * creates a new matrix by appending the vector {@code x} to matrix {@code X} as new row
	 * @param X
	 * @param x
	 * @return
	 */
	public static float[][] appendRows(float[][] X, float[] x){
		float[][] Y = new float[X.length + 1][];
		for (int i = 0; i < X.length; i++) {
			Y[i] = X[i];
		}
		Y[Y.length - 1] = x; 
		return Y;
	}
		
	/**
	 * creates a new matrix by appending the vector {@code x} to matrix {@code X} as new row
	 * @param X
	 * @param x
	 * @return
	 */
	public static long[][] appendRows(long[][] X, long[] x){
		long[][] Y = new long[X.length + 1][];
		for (int i = 0; i < X.length; i++) {
			Y[i] = X[i];
		}
		Y[Y.length - 1] = x; 
		return Y;
	}
	
	/**
	 * creates a new matrix by appending the vector {@code x} to matrix {@code X} as new row
	 * @param X
	 * @param x
	 * @return
	 */
	public static int[][] appendRows(int[][] X, int[] x){
		int[][] Y = new int[X.length + 1][];
		for (int i = 0; i < X.length; i++) {
			Y[i] = X[i];
		}
		Y[Y.length - 1] = x; 
		return Y;
	}
	
	/**
	 * creates a new matrix by appending the vector {@code x} to matrix {@code X} as new row
	 * @param X
	 * @param x
	 * @return
	 */
	public static boolean[][] appendRows(boolean[][] X, boolean[] x){
		boolean[][] Y = new boolean[X.length + 1][];
		for (int i = 0; i < X.length; i++) {
			Y[i] = X[i];
		}
		Y[Y.length - 1] = x; 
		return Y;
	}
	

	/**
	 * creates a new matrix by appending the vector {@code x} to matrix {@code X} as new row
	 * @param X
	 * @param x
	 * @return
	 */
	public static byte[][] appendRows(byte[][] X, byte[] x){
		byte[][] Y = new byte[X.length + 1][];
		for (int i = 0; i < X.length; i++) {
			Y[i] = X[i];
		}
		Y[Y.length - 1] = x; 
		return Y;
	}
	
	/**
	 * creates a new matrix by appending the vector {@code x} to matrix {@code X} as new row
	 * @param X
	 * @param x
	 * @return
	 */
	public static Object[][] appendRows(Object[][] X, Object[] x){
		Object[][] Y = new Object[X.length + 1][];
		for (int i = 0; i < X.length; i++) {
			Y[i] = X[i];
		}
		Y[Y.length - 1] = x; 
		return Y;
	}
	
	/**
	 * creates a new matrix by appending the vector {@code x} to matrix {@code X} as new row
	 * @param X
	 * @param x
	 * @return
	 */
	public static String[][] appendRows(String[][] X, String[] x){
		String[][] Y = new String[X.length + 1][];
		for (int i = 0; i < X.length; i++) {
			Y[i] = X[i];
		}
		Y[Y.length - 1] = x; 
		return Y;
	}
	
	/**
	 * appends the matrix {@code X} by matrix {@code Y} as new columns
	 * <br>if dimensions match
	 * @param X
	 * @param Y
	 * @return
	 */
	public static double[][] appendColumns(double[][] X, double[][] Y){
		int c1 = X[0].length;
		int r1 = X.length;
		int c2 = Y[0].length;
		int r2 = Y.length;
		if (r1 == r2) {
			double[][] Z = new double[r1][c1 + c2];
			for (int n = 0; n < r1; n++) {
				System.arraycopy(X[n], 0, Z[n], 0, c1);
				System.arraycopy(Y[n], 0, Z[n], c1, c2);
			}
			return Z;
		} else {
			throw new IllegalArgumentException("matrices X" + matrixDimensionsToString(X) + " and Y" + matrixDimensionsToString(Y) + " do not have matching columns for appending by row!");
		}
	}
	
	
	/**
	 * enforces symmetric values of the upper triangular of {@code X} to the rest of the matrix
	 * <br>only works with quadratic matrices
	 * <br>the method changes the matrix itself, it does not create a copy
	 * @param X
	 * @return
	 * @throws IllegalArgumentException when the matrix {@code X} is not square
	 */
	public static void makeSymmetric(double[][] X){
		if (!isSquare(X)) {
			throw new IllegalArgumentException("The input matrix must be squared");
		}
		
		int rows = X.length;
		int cols = X[0].length;
		
		for (int i = 0; i < rows - 1; i++) {
			for (int j = 0 + i; j < cols; j++) {
				X[j][i] = X[i][j];
			}		
		}		
	}
		
	public static String matrixDimensionsToString(double[][] X) {
		StringBuilder sb = new StringBuilder();
		int n1 = X[0].length;
		int m1 = X.length;
		sb.append("[").append(m1).append("]").append("[").append(n1).append("]");
		return sb.toString();
	}
	
	public static String matrixDimensionsToString(char[][] X) {
		StringBuilder sb = new StringBuilder();
		int n1 = X[0].length;
		int m1 = X.length;
		sb.append("[").append(m1).append("]").append("[").append(n1).append("]");
		return sb.toString();
	}
	
	public static String matrixDimensionsToString(boolean[][] X) {
		StringBuilder sb = new StringBuilder();
		int n1 = X[0].length;
		int m1 = X.length;
		sb.append("[").append(m1).append("]").append("[").append(n1).append("]");
		return sb.toString();
	}

	public static String matrixDimensionsToString(String[][] X) {
		StringBuilder sb = new StringBuilder();
		int n1 = X[0].length;
		int m1 = X.length;
		sb.append("[").append(m1).append("]").append("[").append(n1).append("]");
		return sb.toString();
	}
		
	public static String matrixDimensionsToString(byte[][] X) {
		StringBuilder sb = new StringBuilder();
		int n1 = X[0].length;
		int m1 = X.length;
		sb.append("[").append(m1).append("]").append("[").append(n1).append("]");
		return sb.toString();
	}
	
	public static String matrixDimensionsToString(float[][] X) {
		StringBuilder sb = new StringBuilder();
		int n1 = X[0].length;
		int m1 = X.length;
		sb.append("[").append(m1).append("]").append("[").append(n1).append("]");
		return sb.toString();
	}
	
	public static String matrixDimensionsToString(long[][] X) {
		StringBuilder sb = new StringBuilder();
		int n1 = X[0].length;
		int m1 = X.length;
		sb.append("[").append(m1).append("]").append("[").append(n1).append("]");
		return sb.toString();
	}
	
	public static String matrixDimensionsToString(int[][] X) {
		StringBuilder sb = new StringBuilder();
		int n1 = X[0].length;
		int m1 = X.length;
		sb.append("[").append(m1).append("]").append("[").append(n1).append("]");
		return sb.toString();
	}
			
	private static  void checkMatrixProdDimensions(double[][] X, double[][] Y) {
		if (X[0].length != Y.length) {
            throw new IllegalArgumentException("Invalid matrix dimensions for multiplication, " + matrixDimensionsToString(X) + " --> " + matrixDimensionsToString(Y));
		}
	}
	
	private static  void checkMatrixProdDimensions(double[][] A, double[] x) {
		if (A[0].length != x.length) {
            throw new IllegalArgumentException("Invalid dimensions for multiplication, " + matrixDimensionsToString(A) + " --> " + x.length);
		}
	}
	
}
