package net.sytes.botg.array.math;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.IntStream;

public class Parallel {

	private static final int NUM_THREADS = Runtime.getRuntime().availableProcessors();
    private static final int TILE_SIZE = 32; // Adjust the tile size based on your system and matrix sizes
	    
    /**
	 * multiplies two double[][] matrices in parallel using ExecutorService and Runnable
	 * @param matrix1
	 * @param matrix2
	 * @return
	 * @throws Exception
	 */
    public static double[][] product(double[][] matrix1, double[][] matrix2) throws Exception {
        int rows1 = matrix1.length;
        int cols1 = matrix1[0].length;
        int rows2 = matrix2.length;
        int cols2 = matrix2[0].length;

        if (cols1 != rows2) {
            throw new IllegalArgumentException("Incompatible matrix dimensions");
        }

        double[][] result = new double[rows1][cols2];

        ExecutorService executorService = Executors.newFixedThreadPool(NUM_THREADS);

        try {
            Future<?>[] futures = new Future[rows1];

            for (int i = 0; i < rows1; i++) {
                futures[i] = executorService.submit(new MatrixRowMultiplier(matrix1, matrix2, result, i));
            }

            for (int i = 0; i < rows1; i++) {
                futures[i].get(); // Wait for each row multiplication to complete
            }

        } finally {
            executorService.shutdown();
        }

        return result;
    }
    
	/**
	 * multiplies two int[][] matrices in parallel
	 * @param matrix1
	 * @param matrix2
	 * @return
	 * @throws Exception
	 */
    public static int[][] multiplyMatrices(int[][] matrix1, int[][] matrix2) throws Exception {
        int rows1 = matrix1.length;
        int cols1 = matrix1[0].length;
        int rows2 = matrix2.length;
        int cols2 = matrix2[0].length;

        if (cols1 != rows2) {
            throw new IllegalArgumentException("Incompatible matrix dimensions");
        }

        int[][] result = new int[rows1][cols2];

        ExecutorService executorService = Executors.newFixedThreadPool(NUM_THREADS);

        try {
            Future<?>[] futures = new Future[rows1];

            for (int i = 0; i < rows1; i++) {
                futures[i] = executorService.submit(new IntMatrixRowMultiplier(matrix1, matrix2, result, i));
            }

            for (int i = 0; i < rows1; i++) {
                futures[i].get(); // Wait for each row multiplication to complete
            }

        } finally {
            executorService.shutdown();
        }

        return result;
    }
    
	/**
	 * Parallel Matrix multiplication using Streams
	 * @param X first matrix 'm×n'
	 * @param Y second matrix 'n×p'
	 * @return result matrix 'm×p'
	 */
	public static double[][] product2(double[][] X, double[][] Y) {
	    int m = X.length; // rows of 'a' matrix
	    int n = X[0].length; // columns of 'a' matrix and rows of 'b' matrix
	    int p = Y[0].length; // columns of 'b' matrix
		return IntStream.range(0, m)
	            .parallel() // comment this line to check the sequential stream
	            .mapToObj(i -> IntStream.range(0, p)
	                    .mapToDouble(j -> IntStream.range(0, n)
	                            .mapToDouble(k -> X[i][k] * Y[k][j])
	                            .sum())
	                    .toArray())
	            .toArray(double[][]::new);
	}

    /**
     * is used in {@code multiplyMatrices} as Runnable Task to be run in parallel threads
     * @author hillenbrand
     *
     */
    private static class IntMatrixRowMultiplier implements Runnable {
        private final int[][] matrix1;
        private final int[][] matrix2;
        private final int[][] result;
        private final int row;

        public IntMatrixRowMultiplier(int[][] matrix1, int[][] matrix2, int[][] result, int row) {
            this.matrix1 = matrix1;
            this.matrix2 = matrix2;
            this.result = result;
            this.row = row;
        }

        @Override
        public void run() {
            int cols1 = matrix1[0].length;
            int cols2 = matrix2[0].length;

            for (int kStart = 0; kStart < cols1; kStart += TILE_SIZE) {
                int kEnd = Math.min(kStart + TILE_SIZE, cols1);

                for (int jStart = 0; jStart < cols2; jStart += TILE_SIZE) {
                    int jEnd = Math.min(jStart + TILE_SIZE, cols2);

                    for (int k = kStart; k < kEnd; k++) {
                        for (int j = jStart; j < jEnd; j++) {
                            result[row][j] += matrix1[row][k] * matrix2[k][j];
                        }
                    }
                }
            }
        }
    }

    /**
     * is used in {@code multiplyMatrices} as Runnable Task to be run in parallel threads
     * @author hillenbrand
     *
     */
    private static class MatrixRowMultiplier implements Runnable {
        private final double[][] matrix1;
        private final double[][] matrix2;
        private final double[][] result;
        private final int row;

        public MatrixRowMultiplier(double[][] matrix1, double[][] matrix2, double[][] result, int row) {
            this.matrix1 = matrix1;
            this.matrix2 = matrix2;
            this.result = result;
            this.row = row;
        }

        @Override
        public void run() {
            int cols1 = matrix1[0].length;
            int cols2 = matrix2[0].length;

            for (int kStart = 0; kStart < cols1; kStart += TILE_SIZE) {
                int kEnd = Math.min(kStart + TILE_SIZE, cols1);

                for (int jStart = 0; jStart < cols2; jStart += TILE_SIZE) {
                    int jEnd = Math.min(jStart + TILE_SIZE, cols2);

                    for (int k = kStart; k < kEnd; k++) {
                        for (int j = jStart; j < jEnd; j++) {
                            result[row][j] += matrix1[row][k] * matrix2[k][j];
                        }
                    }
                }
            }
        }
    }
    
}
