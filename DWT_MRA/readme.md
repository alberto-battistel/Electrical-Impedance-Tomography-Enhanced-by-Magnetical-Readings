

# Interpolation of Discrete Wavelet Basis and Multiresolution Basis on a FEM Mesh

The winner is **Haar discrete wavelet**

## `MRA_polar_mesh_test.m`

**2D**

- Plots the Haar multiresolution wavelet.
- Interpolates the wavelets on the FEM mesh using polar coordinates.
- Interpolates a fictional distribution with the wavelets using different approaches.
- **Best result:** **Lasso**, as it is numerically stable and gives a small \( L_2 \) norm.
- Basis size = 256, (n*n)

### Approaches

#### Least Squares
Not numerically stable
- \( L_2 = 4.975434 \) with \( n = 16 \)

#### Inner Product
- \( L_2 = 36.855245 \) with \( n = 16 \)

#### Proper Orthogonal Reduction?
- \( L_2 = 5.602834 \) with \( n = 16 \)
- \( l = 232 \)

#### Proper Orthogonal Reduction with New Basis
- \( L_2 = 7.355316 \) with \( n = 16 \)

#### Lasso
- \( L_2 = 5.014755 \) with \( n = 16 \)

#### Lasso with Proper Reduction
- \( L_2 = 9.312381 \) with \( n = 16 \)





## `dwt_polar_mesh_test.m`

**2D**

- Plots the Haar discrete wavelet.
- Interpolates the wavelets on the FEM mesh using polar coordinates.
- Interpolates a fictional distribution with the wavelets using different approaches.
- **It works very well.**
- **Basis size:** 256, ~(2n * 2n)

### Approaches

#### Least Squares
- \( L_2 = 4.975434 \) with \( n = 8 \)

#### Inner Product
- \( L_2 = 36.855245 \) with \( n = 8 \)

#### Lasso
- \( L_2 = 4.975514 \) with \( n = 8 \)


## `dwt_cylindrical_mesh_test.m`

- Plots the Haar discrete wavelet.
- Interpolates the wavelets on the FEM mesh using cylindrical coordinates.
- Interpolates a fictional distribution with the wavelets using different approaches.
- **It works very well.**
- **Basis size:** 4096, ~(2n1 * 2n2 * 2n3)

### Approaches

#### Least Squares
It is numerically stable, cond = 587.4037
- \( L_2 = 43.290449 \) with \( n = 4096 \)

#### Lasso
- \( L_2 = 43.290729 \) with \( n = 4096 \)



## `dwt_polar_mesh_test2.m`

**2D**

Similar to `dwt_polar_mesh_test.m` but different mapping

- Plots the Haar discrete wavelet.
- Interpolates the wavelets on the FEM mesh using polar coordinates, but instead of mapping on r and theta it uses 2r and theta/2.
- Interpolates a fictional distribution with the wavelets using different approaches.
- **Only lasso works well**
- **Basis size:** 256, ~(2n * 2n)
- sure problems with mapping

### Approaches

#### Least Squares
numerically unstable
- \( L_2 = 5.087603 \) with \( n = 8 \)

#### Inner Product
- \( L_2 = 35.120789 \) with \( n = 8 \)

#### Lasso
- \( L_2 = 5.098968 \) with \( n = 8 \)




## `dwt_test.m`

**1D**

plot and interpolate a function in 1D with Haar discrete wavelets.


## `dwt_mesh_test`

**2D**

- Interpolate Haar wavelets on mesh from a rectangular plane.
- Interpolates a fictional distribution with the wavelets using different approaches.

It works but a lot of the wavelets are wasted on empty space.
