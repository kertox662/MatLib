# MatLib
A simple C++ Matrix Library. This was made to learn about Linear Algebra.

# What Does it Do?
The library currently contains several classes:
 - Matrix<int,int>
 - Vec3
 - Vec4
 
 These classes can be used as storage classes, as well as be transformed. The library provides functions to create different transformation matrices.
 
# On the Matrix Class
 - To construct a Matrix, a default constructor can be used to fill the matrix with 0s, or a 2D-array of the same size can be provided as a parameter.
 - Accessing elements in the Matrix is the same as a 2D-array.
 - To add or multtiply a matrix the `add(Matrix<n,m>)` and `mult(Matrix<n,m>)` member functions can be called, passing in a second the Matrix. Another way to do this is to use the `+` or `*` operators.
   - NOTE: The responsibility for proper multiplication and addition is on the user. The library will not allow for undefined behaviour during compilation.
 - For square matrices, certain other functions are defined.
    - The static member `Identity()` returns an Indentity matrix of the same size.
    - `determinant()` returns the value of the determinant.
    - `inverse()` and `inverse(int*)` gets the inverse of the matrix if it exists. The second version of the member puts an error code into the pointer. `0` means success and `1` means that the matrix is degenerate.
    - Certain helper functions that are used in `inverse()` like `getCofactors()` and `getAdjugate()` which return the type of matrices of the same dimension.

# On the Vector Classes
 - Vec3 can be made by passing in 3 `double` values or a `double[3]`.
 - It provides basic vector operations:
    - Addition in `add(Vec3)` or the `+` operator.
    - Dot product in `dot(Vec3)`
    - Cross product in `cross(Vec3)`.
    - Magnitude in `mag()`
    - The angle between vectors in `angleBetween(Vec3)`
 - A Vec3 can be transformed by a 3x3 matrix by calling `transform(Matrix<3,3>)` or a 4x4 by calling `transform(Matrix<4,4>)`.
 - NOTE: Although Vec4 can be used as a 4D vector, that is not the intended purpose in this library. It is there to allow for translation by a translation matrix.
 
# On Transforming 
The libary provides a namespace `Transform` where certain aspects of transformation are available to the user.

- `make2DRot(double)` makes a 2x2 matrix with that rotates by the angle provided (in radians).
- `make[Axis]Rot(double)` makes a 3x3 matrix. `[Axis]` must be replaced by either X,Y or Z. The resulting matrix rotates about the axis by the angle provided.
- `makeRot(double, double, double)` gives the product of the 3 rotation matrices.
- `makeTranslation(Vec3)` gives a 4x4 matrix that translates a vector by the given Vec3.
- `makeScale(double)` gives a 3x3 scaling matrix that scales by the given double.
