# Unit-Tests for Load-Balanced Chemistry Model

This folder contains tests and test cases for unit tests of the load-balanced 
chemistry models.

Test cases can be found in `Cases` and the source code
for the tests in the `src` folder.

## Compilation and Usage

The tests use Catch2 as a driving library, included as
a git sub-module.

### Prerequisites

1. CMake installed (Minimum version 3.5)
2. Update git-submodule with:
   `git submodule init && git submodule update`

### Compile and Run Tests

1. Execute `./Allwmake` in the tests directory.
   This will compile and install Catch2 locally in your
   repository and compile all tests with `wmake`
2. Run `./Allrun` to run all tests
3. Alternative: Go to the test case and execute the respective
   test with the correct tag (For tags see Catch2 manual).
   An example for testing the chemistry models in parallel is:
   ```bash
cd Cases/Case-chemistry
blockMesh
decomposePar 
mpirun -np 4 ../../unitTests.exe [chemistry] --parallel
   ```

