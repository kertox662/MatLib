# MatLib
A simple C++ Matrix Library. This was made to learn about Linear Algebra.

## What Does it Do?
The library currently contains several classes:
 - Matrix<int,int>
 - vec3
 - vec4
 
 These classes can be used as storage classes, as well as be transformed. The library provides functions to create different transformation matrices. **These classes, as well as everything in this library is in the _Matlib_ namespace.**
 
## On the Matrix Class
 - To construct a Matrix, a default constructor can be used to fill the matrix with 0s, or a 2D-array of the same size can be provided as a parameter.
 - Accessing elements in the Matrix is the same as a 2D-array.
 - To add or multtiply a matrix the `add(Matrix<n,m>)` and `mult(Matrix<n,m>)` member functions can be called, passing in a second the Matrix. Another way to do this is to use the `+` or `*` operators.
   - NOTE: The responsibility for proper multiplication and addition is on the user. The library will not allow for undefined behaviour during compilation.
 - For square matrices, certain other functions are defined.
    - The static member `Identity()` returns an Indentity matrix of the same size.
    - `determinant()` returns the value of the determinant.
    - `inverse()` and `inverse(int*)` gets the inverse of the matrix if it exists. The second version of the member puts an error code into the pointer. `0` means success and `1` means that the matrix is degenerate.
    - Certain helper functions that are used in `inverse()` like `getCofactors()` and `getAdjugate()` which return the type of matrices of the same dimension.

## On the Vector Classes
#### vec3
 - vec3 can be made by passing in 3 `double` values or a `double[3]`.
 - An additional way to make vec3 is by passing in a `Matrix<3,1>`.
 - It provides basic vector operations:
    - Addition in `add(vec3)` or the `+` operator.
    - Subtraction in `sub(vec3)` or the `-` operator.
    - Get the opposite vector using `neg()` or the `-` operator.
    - Scalar multiplication in `mult(double)` or the `*` operator.
    - Dot product in `dot(vec3)`
    - Cross product in `cross(vec3)`.
    - Magnitude in `mag()`
    - The angle between vectors in `angleBetween(vec3)`
 - A vec3 can be transformed by a 3x3 matrix by calling `transform(Matrix<3,3>)` or a 4x4 by calling `transform(Matrix<4,4>)`.
 #### vec4
 - vec4 can be made in the same way as vec3, with the addition of using vec3. In that case, the fourth value will be set to 1.
 - The following operations are defined:
   - Addition
   - Substraction
   - Opposite Vector
   - Scalar multiplication
   - Dot Product
   - Magnitude
 
## On Transforming 
The libary provides a namespace `Transform` where certain aspects of transformation are available to the user.

- `make2DRot(double)` makes a 2x2 matrix with that rotates by the angle provided (in radians).
- `make[Axis]Rot(double)` makes a 3x3 matrix. `[Axis]` must be replaced by either X,Y or Z. The resulting matrix rotates about the axis by the angle provided.
- `makeRot(double, double, double)` gives the product of the 3 rotation matrices.
- `makeTranslation(vec3)` gives a 4x4 matrix that translates a vector by the given vec3.
- `makeScale(double)` gives a 3x3 scaling matrix that scales by the given double.
