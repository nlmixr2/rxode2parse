# Cran comments
- Fixes the m1mac issues flagged by CRAN, including the overlooked issue on 1323:19
- Removes solved systems with gradients for intel c++ compilers (because they didn't compile)
  - Have updated reverse dependencies to check if linear compartment solutions are present before running tests
