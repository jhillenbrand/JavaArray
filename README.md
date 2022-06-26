# JavaArray

Java project containing classes with static methods to manipulate arrays, doing math operations,  signal processing and more ...

## Packages

### array

* [ArUtils](src/net/sytes/botg/array/ArUtils.java)
* [ConvertArray](src/net/sytes/botg/array/ConvertArray.java)
* [SortArray](src/net/sytes/botg/array/SortArray.java)
* [SearchArray](src/net/sytes/botg/array/SearchArray.java)

### math

and Math functions for signal processing and feature extractions and many more ...
* [Scalar](src/net/sytes/botg/array/math/Scalar.java)
* [Vec](src/net/sytes/botg/array/math/Vec.java)
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

The repo also contains code samples for MATLAB and C (compiled under Windows10), that compare execution time with the UnitTest_ScalarProd:
<br>for comparison the scalar product of two double arrays with 10e6 elements is computed
* Java		-> 20.000 - 35.000 �s
* C    		-> 25.000 - 45.000 �s
* MATLAB	-> 300.000 - 1e6 �s

<br>I'm not claiming that my Java code is faster than native C, because probably the implementation is not using all advantages of native C, but it shows me that there's hope for me doing my Signal Processing in Java (:P)

## Releases
this repo contains *.jar releases for every release/snapshot branch (except master), there will be a downloadable *.jar, that can be intergrated into your project, if you do not wish to compile the src yourself.

## Licenses used
this Project: MIT

code used within Project is subject to EPL 2.0
