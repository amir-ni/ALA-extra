# Finding the orthogonal projection of a point onto an affine subspace

This repository contains an implementation of an iterative method to find the orthogonal projection of a given point to the solution set of a system of linear equations in `Julia`. This method can also be used to solve a system of linear equations.

This implementation is done as an extra credit project for the [Applied Linear Algebra course](http://sharif.edu/~mtefagh/ala/home.html).

## Usage

The method is available as a function in `src/solver.jl`.

## Reference

The implementation is based on the method proposed in [this paper](https://www.sciencedirect.com/science/article/pii/S0024379506005040).

## Experiments

In the reference paper, the implementation was done using MATLAB, and some numerical results were presented. To compare the performance, numerical experiments were repeated using this implementation, and the result are available below. Results obtained on a laptop with Intel i7-7700HQ @ 2.80GHz processor. These results are reproducible using `src/experiments.jl`.

We also considered implementation remarks made in the paper and implemented them in our code.

### Method

The below description is based on the evaluation method used in the original paper.

We used random matrixes and Hilbert matrixes to evaluate the method. In the former case, we shifted a random matrix to get also negative elements: A = rand(m, n) − 0.5. In the latter case we generated a Hilbert matrix A with A(i, j) = 1/(i + j − 1). In all cases, we took a vector z of 1s to be a solution and defined b = Az. The presented results concern the ∞-norm.

The below tables are equivalent to tables 2,3,4, and 5 in the reference paper, respectively.

### Finding projections: numerical results for random matrices

m | n | ∥y-x₀∥ | ∥b-Ay∥/∥b∥ | Time
--- | --- | --- | --- |---:
1000 | 1000 | 1.0e+0 | 4.1e-12 | 25.27
500 | 500 | 1.0e+0 | 6.6e-11 | 1.47
200 | 200 | 1.0e+0 | 1.0e-12 | 0.09
100 | 100 | 1.0e+0 | 2.4e-13 | 0.04
1000 | 2000 | 2.0e+0 | 4.6e-13 | 51.30
500 | 2000 | 1.8e+0 | 3.1e-15 | 12.37
200 | 2000 | 1.2e+0 | 1.1e-15 | 1.17
100 | 2000 | 7.4e-1 | 9.9e-16 | 0.28
50 | 2000 | 5.3e-1 | 1.4e-15 | 0.07
20 | 2000 | 2.4e-1 | 4.2e-16 | 0.01
10 | 2000 | 2.5e-1 | 4.9e-16 | 0.00
200 | 10000 | 5.2e-1 | 3.4e-15 | 9.88
200 | 20000 | 4.4e-1 | 2.9e-15 | 37.82
200 | 50000 | 2.8e-1 | 6.6e-15 | 102.42
100 | 20000 | 3.2e-1 | 2.5e-15 | 7.12
100 | 100000 | 1.5e-1 | 9.1e-15 | 49.24
100 | 200000 | 9.1e-2 | 6.7e-15 | 104.54
50 | 100000 | 9.0e-2 | 7.8e-15 | 11.10
50 | 200000 | 7.0e-2 | 2.2e-14 | 26.96
50 | 400000 | 6.1e-2 | 8.9e-15 | 50.53
20 | 200000 | 4.3e-2 | 9.0e-16 | 3.92
20 | 400000 | 3.1e-2 | 1.1e-15 | 7.34
20 | 1000000 | 1.5e-2 | 2.e-15 | 21.18
10 | 100000 | 2.1e-2 | 1.0e-15 | 0.45
10 | 200000 | 2.3e-2 | 8.3e-16 | 1.01
10 | 400000 | 1.7e-2 | 7.3e-16 | 1.76
10 | 1000000 | 1.4e-2 | 2.9e-15 | 4.42
10 | 2000000 | 8.7e-3 | 3.0e-15 | 9.38

### Finding projections: numerical results for Hilbert matrices

m | n | ∥y-x₀∥ | ∥b-Ay∥/∥b∥ | Time
--- | --- | --- | --- |---:
200 | 200 | 1.98e+0 | 3.52e-7 | 0.81
100 | 100 | 1.00e+0 | 6.61e-10 | 0.03
50 | 2000 | 2.52e+1 | 6.02e-6 | 0.14
10 | 2000 | 4.75e+0 | 2.01e-8 | 0.00
200 | 10000 | 5.84e+1 | 3.76e-5 | 12.77
200 | 20000 | 2.22e+2 | 9.41e-5 | 36.72
100 | 20000 | 5.37e+2 | 2.22e-5 | 7.46
100 | 100000 | 5.60e+2 | 4.35e-5 | 45.33
50 | 100000 | 2.90e+1 | 8.57e-6 | 9.77
50 | 200000 | 3.12e+1 | 5.28e-7 | 23.93
20 | 200000 | 2.06e+1 | 3.34e-7 | 3.45
20 | 400000 | 1.86e+1 | 3.00e-7 | 6.66
10 | 100000 | 1.08e+1 | 1.02e-7 | 0.41
10 | 200000 | 1.14e+1 | 6.56e-7 | 0.87
10 | 400000 | 2.28e+1 | 4.35e-7 | 1.75
10 | 1000000 | 1.64e+1 | 1.91e-7 | 4.20
10 | 2000000 | 2.00e+1 | 3.66e-7 | 7.86

### Linear equations: numerical results for random matrices

m | n | ∥x-z∥ | ∥b-Ax∥/∥b∥ | Time
--- | --- | --- | --- |---:
100 | 100 | 1.00e+0 | 2.25e-13 | 0.60
200 | 200 | 1.00e+0 | 1.16e-12 | 0.23
1000 | 1000 | 1.00e+0 | 6.60e-13 | 22.02
2000 | 2000 | 1.00e+0 | 6.85e-11 | 210.07
100 | 100000 | 1.40e-1 | 6.10e-15 | 39.53
10 | 1000000 | 9.04e-3 | 3.92e-15 | 4.71

### Linear equations: numerical results for Hilbert matrices

m | n | ∥x-z∥ | ∥b-Ax∥/∥b∥ | Time
--- | --- | --- | --- |---:
10 | 10 | 1.00e+0 | 1.13e-12 | 0.00
20 | 20 | 1.00e+0 | 3.98e-11 | 0.00
50 | 50 | 1.00e+0 | 8.67e-11 | 0.00
100 | 100 | 1.00e+0 | 6.61e-10 | 0.09
200 | 200 | 1.98e+0 | 3.52e-7 | 0.37
100 | 100000 | 5.60e+2 | 4.35e-5 | 39.73
10 | 1000000 | 1.64e+1 | 1.91e-7 | 4.32

## License

Distributed under the MIT License. See `LICENSE` for more information.

## Contact

Amirhossein Nadiri - amirhossein@nadiri.me

Project Link: [https://github.com/amir-ni/ALA-extra](https://github.com/amir-ni/ALA-extra)