Optimet
=======

[![homepage](https://img.shields.io/badge/homepage-url-blue.svg)](https://www.ee.ucl.ac.uk/~npanoiu/Software.html)
[![license](https://img.shields.io/badge/license-GPL-blue.svg)](LICENSE)
[![manual](https://img.shields.io/badge/manual-docx-yellow.svg)](manuals/manual.docx)
[![manual](https://img.shields.io/badge/manual-pdf-yellow.svg)](manuals/manual.pdf)

Optimet can be used for the scattering analysis of light at fundamental and second-harmonic frequencies on distribution of homogeneous spherical or nonspherical
nanoparticles embedded in a homogeneous medium. It can run in parallel on large clusters, using MPI.
It accepts several linear-system solvers, including a direct solver from
[Scalapack](http://www.netlib.org/scalapack/), and iterative solvers from
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

- [CMake](https://cmake.org/): The build system. Must be installed independently.
- MPI: Required to run in parallel. Must be installed independently. Essential for nonspherical particles.
- Scalapack: (optional) Parallel linear algebra. Must be installed independently. Only useful when
  compiling with MPI.
- [Belos](https://trilinos.org/packages/belos/): (optional) A library of iterative solvers. Must be
  installed independently. Only useful when compiling with MPI. Belos comes as a part of Trilinos package. Version 12-10-1 of Trilinos is used (git checkout trilinos-release-12-10-1).
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

Installation
------------

The build system is [CMake](https://cmake.org/). It can generate a build environment from a fair
number of systems, from Unix makefiles to Xcode projects.

The canonical usage is as follows:

```
cd /path/to/optimet_source
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -Ddompi=ON -Ddoarshp=OFF ..
make
```

The example above will *not* compile the code responsible for the analysis of nonspherical targets. Set this option
to `ON` if you would rather have it. The code used for the TMM analysis of nonspherical particles has to be compiled with the MPI option set to ON.

The executable `Optimet3D` should be directly in the build directory.

Supported Platforms
-------------------

Optimet is routinely tested on:

- Linux
- MacOS

It is known to compile with:

- gnu-g++ > 4.8
- clang
- intel > 13.0
