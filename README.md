# meshing-circle
solving Two-Dimensional Poisson Equations with Dirichlet boundary, ![1](https://latex.codecogs.com/gif.latex?%5COmega%20%3D%20%5C%7B%20%28x%2Cy%29%20%7Cx%5E%7B2%7D&plus;y%5E%7B2%7D%20%5Cleq%20R%5C%7D)
# Using Dinite Element Methods (FEM)

cmesh.py - create and interact with the mesh

FEM - sovle for ![2](https://latex.codecogs.com/gif.latex?-%5CDelta%20u%3D%20f)

Fem_2 - sovle for general form: 

![3](https://latex.codecogs.com/gif.latex?%5Cdpi%7B120%7D%20%5Clarge%20-%5Ba_%7B00%7D%28x%2C%20y%29%20%5Cfrac%7B%5Cpartial%5E%7B2%7D%20u%7D%7B%7B%5Cpartial%20x%7D%5E%7B2%7D%7D%20&plus;%20a_%7B01%7D%28x%2C%20y%29%20%5Cfrac%7B%5Cpartial%5E%7B2%7D%20u%7D%7B%5Cpartial%20x%20%7B%5Cpartial%20y%7D%7D%20&plus;%20a_%7B10%7D%28x%2C%20y%29%20%5Cfrac%7B%5Cpartial%5E%7B2%7D%20u%7D%7B%5Cpartial%20y%20%7B%5Cpartial%20x%7D%7D%20&plus;%20a_%7B11%7D%28x%2C%20y%29%20%5Cfrac%7B%5Cpartial%5E%7B2%7D%20u%7D%7B%7B%5Cpartial%20y%7D%5E%7B2%7D%7D%20%5D%20&plus;%20a_%7B0%7D%28x%2Cy%29.u%20%3D%20f%28x%2Cy%29)
