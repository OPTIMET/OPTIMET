Optimet
=======

[![homepage](https://img.shields.io/badge/homepage-url-blue.svg)](https://www.ee.ucl.ac.uk/~npanoiu/Software.html)
[![license](https://img.shields.io/badge/license-GPL-blue.svg)](LICENSE)
[![manual](https://img.shields.io/badge/manual-docx-yellow.svg)](manuals/manual.docx)
[![manual](https://img.shields.io/badge/manual-pdf-yellow.svg)](manuals/manual.pdf)

Optimet is a simulation of multiple-scattering of light on a distribution of homogeneous spherical
nanoparticles embedded in a homogeneous medium. It can run in parallel on large clusters, using MPI.
It accepts several linear-system solvers, including a direct solver from
[Scalapack](http://www.netlib.org/scalapack/), and the iterative solvers from
[Belos](https://trilinos.org/packages/belos/). Optimet is a c++11 program.

Optimet is distributed under the GNU Public License.

Usage
=====

Example inputs are provided in the `example` subdirectory. The input is in XML format. It can
specify any arbitrary distribution of particles.

The compiled program is called with a single argument, the path to the XML input:

```
Optimet3D input.xml
```

A manual is available in [docx](manuals/manual.docx) and [pdf](manuals/manual.pdf) formats.

Installation
============

Dependencies
------------

Optimet depends on the following external packages. In most cases, the build system will try to
download and install dependencies it cannot find on the system.

- [CMake](https://cmake.org/): The build system. Must be installed independantly.
- MPI: Required to run in parallel. Must be installed independantly.
- Scalapack: (optional) Parallel linear algebra. Must be installed independantly. Only usefull whhen
  compiling with MPI.
- [Belos](https://trilinos.org/packages/belos/): (optional) A library of iterative solvers. Must be
  installed independantly. Only usefull when compiling with MPI.
- [Boost](http://www.boost.org/): (required) A set of peer-reviewed c++ libraries. At this juncture,
  Optimet only requires the [math special functions
  module](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/special.html). Automatically
  downloaded if unavailable.
- [Eigen](http://eigen.tuxfamily.org/Findex.php?title=Main_Page): (required) A c++ linear algebra
  library.  Automatically downloaded if unavailable.
- [hdf5](https://support.hdfgroup.org/HDF5/): (required) A standard library for handling data.
  Automatically downloaded if unavailable. Only the C bindings are used. The c++ bindings are not
  necessary.
- [GSL](https://www.gnu.org/software/gsl/): (required) Numerical library for C and C++.
  Automatically downloaded if unavailable.
- [F2C](http://www.netlib.org/f2c/): (required) Fortran to C library. Required by the Amos package.
- [Benchmark](https://github.com/google/benchmark): (optional) Micro-benchmarking library from
  google. Automatically downloaded if unavailable *and* benchmarks are requested.
- [Catch](https://github.com/philsquared/Catch): (optional) A c++ testing framework. Downloaded
  automatically if tests are required.

Installation
------------

The build system is [CMake](https://cmake.org/). It can generate a build environment from a fair
number of systems, from Unix makefiles to Xcode projects.

The cannonical usage is as follows:

```
cd /path/to/optimet_source
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -Ddompi=ON -Ddotesting=OFF -Ddobenchmarks=OFF ..
make
```

The example above will *not* compile the tests and benchmarks. Leave these options out (or set them
to `ON`) if you would rather have them. Similarly, the example compiles with MPI parallelization
enabled. Remove or set to `ON` to compile a serial code only.

The executable `Optimet3D` should be directly in the build directory. We currently provide not
installation mechanism.

The benchmarks and unit-tests (when compiled) are in the `benchmarks` and `unittests` subdirectories
of the build directory. They consists in a number of executables which can be invoked manually.
Alternatively, the tests can be invoked from the build directory with `make test` or `ctest .`.

Supported Platforms
-------------------

Optimet is routinely tested on:

- Linux
- MacOS

It is known to compile with:

- gnu-g++ > 4.8
- clang
- intel > 13.0
