# JavaArray
Java project containing classes with static methods to manipulate arrays, doing math operations, signal processing, geometry and more ...

## Packages

### array

* [ConvertArray](src/net/sytes/botg/array/ConvertArray.java)
  <br><i>wrapping and concatenation for arrays</i>
  ```java
  public static Double[] wrap(final double[] array)
  public static double[] unwrap(final Double[] array)
  public static double[] toDouble(int[] intAr)
  public static double[] shortToDouble(short[] shortAr)
  ...
  ```

### math

and Math functions for signal processing and feature extractions and many more ...
* [Scalar](src/net/sytes/botg/array/math/Scalar.java)
  <br><i>methods that operate on scalar input</i>
  ```java
  public static double roundToDecimals(double x, int decimals)
  public static int randInt(int min, int max)
  public static int closestExponentForBase2(int n)
  public static boolean isDivByPow(int n, int m)
  ...
  ```
* [Vec](src/net/sytes/botg/array/math/Vec.java)
  <i><br>Collection of methods that operates on vectors (mainly double[]) as method input
  <br>it contains methods that transform vectors into scalars, vectors into vectors, vectors into matrices,<br>generates vectors, searches in vectors, permutates vectors and more</i>
  ```java
  public static double max(double[] x)
  public static double rms(double[] x)
  public static double[] slidingMax(double[] x, int w, int s)
  public static double[] lowpass(double[] x, double dt, double f_c)
  public static double[] normalize(double[] x)
  public static double[] zscore(double[] x)
  public static double[] downsample(double[] ar, int m, DownsamplingAlgorithm algorithm)
  public static double[] upsample2(double[] x, int n)
  public static int[] findLocalExtrema(double[] x, int minDistance)
  public static double[] diff(double[] ar)
  public static double[] sub(final double[] ar, int s, int e)
  public static double[] linspace(double start, double end, double step)
  public static double[] nan(int n)
  ...
  ```
* [Mat](src/net/sytes/botg/array/math/Mat.java)

### spectrum
// TODO description missing

### geometry
// TODO description missing

## Testing

unit tests with examples can be found here for the corresponding package:
* [test](test/)

The repo also contains code samples for MATLAB and C (compiled under Win10), that compare execution time with the UnitTest_ScalarProd:
<br>for comparison the scalar product of two double arrays with 10e6 elements is computed
* Java		-> 20.000 - 35.000 µs
* C    		-> 25.000 - 45.000 µs
* MATLAB	-> 300.000 - 1e6 µs

<br>DISCLAIMER: I'm not claiming that my Java code is faster than native C, because probably the implementation is not using all advantages of native C, but it shows me that there's hope for doing my signal processing in Java (:P)

## Licenses used
this Project: MIT

code used within Project is subject to EPL 2.0

## Releases
If you want to include src sode, that does not change anymore, see the list of releases

## Dependencies
see [pom.xml](pom.xml)