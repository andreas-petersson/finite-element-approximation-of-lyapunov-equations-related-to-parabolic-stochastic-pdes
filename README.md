# Finite Element Approximation of Lyapunov Equations Related to Parabolic Stochastic PDEs

This repository contains the MATLAB code implementation accompanying the paper "Finite element approximation of Lyapunov equations related to parabolic stochastic PDEs". The code provides numerical methods for approximating the solution to operator Lyapunov equations associated with linear stochastic partial differential equations (SPDEs) driven by multiplicative noise, and for computing related quadratic functionals.

## Paper Reference

A. Andersson, A. Lang, A. Petersson, and L. Schroer. *Finite element approximation of Lyapunov equations related to parabolic stochastic PDEs*.

## Code Description

The main file `lyapunov.m` contains the following MATLAB functions:

*   **`compute_lyapunov(a,b,g,r,T,Nh,Nk)`**: Computes the matrix representation `L^N_{h,tau}` of the fully discrete approximation to the solution `L(T)` of the operator Lyapunov equation using a finite element method in space and a semi-implicit Euler scheme in time.
*   **`compute_lyapunov_functional(L,u0)`**: Calculates the approximation `Phi^L_{h,tau}(x) = <L^N_{h,tau} I_h x, I_h x>` of the quadratic functional `Phi(x)` using the matrix computed by `compute_lyapunov`.
*   **`mc(a,b,g,r,T,Nh,Nk,N,x0)`**: Approximates the functional `Phi(x)` using a Monte Carlo method based on simulating sample paths of the fully discrete SPDE approximation.
*   **`sample_HSHE_end_norm(...)`**: Helper function for `mc` that generates one sample path of the discretized SPDE and computes the corresponding functional value for that path. Assumes the noise is white in space.
*   **`stiffness_matrix(h)`**: Computes the finite element stiffness matrix for the 1D Laplacian with Dirichlet boundary conditions.
*   **`mass_matrix(h)`**: Computes the finite element mass matrix for linear elements in 1D.
*   **`generate_wiener_matrix(...)`**: Generates a matrix representation for the discretized spatial noise term.
*   **`cholesky_matrix_h_interval()`**: Helper function used in noise generation.

The file also includes a demo section illustrating how to use the primary functions (`compute_lyapunov_functional` and `mc`).

## Usage

1.  Ensure you have MATLAB installed.
2.  Place `lyapunov.m` in the MATLAB path.
3.  Open `lyapunov.m` in MATLAB.
4.  Modify the parameters (coefficients `a, b, g, r`, time `T`, discretization parameters `Nh, Nk`, Monte Carlo samples `N_mc`, initial condition `x0`) in the "Demo Usage" section as needed.
5.  Run the script. The output will show the functional approximations computed via the Lyapunov method and the Monte Carlo method.
