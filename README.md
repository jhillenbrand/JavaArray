# README is not up-to-date
# JavaArray

Java project containing classes with static methods to manipulate arrays, doing math operations, signal processing, geometry and more ...

## Packages

### array

* [ConvertArray](src/net/sytes/botg/array/ConvertArray.java)
  <br><i>wrapping and concatenation for arrays</i>
  ```java

  ```

### math

and Math functions for signal processing and feature extractions and many more ...
* [Scalar](src/net/sytes/botg/array/math/Scalar.java)
* [Vec](src/net/sytes/botg/array/math/Vec.java)
  ```java
  public static double[] highpass(double[] x, double dt, double f_c)
  public static double[] lowpass(double[] x, double dt, double f_c)
  public static double[] normalize(double[] x)
  public static double[] zscore(double[] x)
  public static double[] downsample(double[] ar, int m, DownsamplingAlgorithm algorithm)
  public static double[] upsample2(double[] x, int n)
  public static int[] findLocalExtrema(double[] x, int minDistance)
  public static double[] diff(double[] ar) 
  ```
* [Vec2Mat](src/net/sytes/botg/array/math/Vec2Mat.java)
* [Vec2Scalar](src/net/sytes/botg/array/math/Vec2Scalar.java)

## Testing

unit tests with examples can be found here:
* [ArUtils](test/array/UnitTest_ArUtils.java)
* [ConvertArray](test/array/UnitTest_ConvertArray.java)
* [SearchArray](test/array/UnitTest_SearchArray.java)
* [SortArray](test/array/UnitTest_SortArray.java)
* [Vec2Vec](test/math/UnitTest_Vec2Vec.java)
* [ScalarProd](test/math/UnitTest_ScalarProd.java)

The repo also contains code samples for MATLAB and C (compiled under Win10), that compare execution time with the UnitTest_ScalarProd:
<br>for comparison the scalar product of two double arrays with 10e6 elements is computed
* Java		-> 20.000 - 35.000 µs
* C    		-> 25.000 - 45.000 µs
* MATLAB	-> 300.000 - 1e6 µs

<br>I'm not claiming that my Java code is faster than native C, because probably the implementation is not using all advantages of native C, but it shows me that there's hope for doing my signal processing in Java (:P)

## Documentation
Please download and see the attached [doc](doc/) folder for an overview of all classes and methods

## Releases
this repo contains *.jar releases for every release/snapshot branch (except master), there will be a downloadable *.jar, that can be intergrated into your project, if you do not wish to compile the src yourself.

## Licenses used
this Project: MIT

code used within Project is subject to EPL 2.0
