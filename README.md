# Discrete Calabi Flow

This project is a basic C++ implementation of the discrete Calabi flow algorithm as proposed by [Su et al., 2019](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.13873). The algorithm computes the discrete conformal map of a given triangle mesh, and this implementation specifically focuses on the disk embedding algorithm presented in the paper.

## Dependencies

This project uses `vcpkg` to manage its dependencies. The following libraries are required:
- `Eigen` (3.4.0)
- `OpenMesh` (10.0)