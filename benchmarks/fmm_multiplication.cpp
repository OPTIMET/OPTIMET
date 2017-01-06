#include "Excitation.h"
#include "Geometry.h"
#include "Tools.h"
#include "Types.h"
#include "mpi/FastMatrixMultiply.h"
#if defined(OPTIMET_MPI) && !defined(OPTIMET_JUST_DO_SERIAL)
#include "mpi/Session.h"
#endif
#include <chrono>
#include <iostream>

namespace {
constexpr optimet::t_real default_wavelength() { return 750e-9; }
constexpr optimet::t_real default_length() { return 2000e-9; }
optimet::Matrix<optimet::t_real> fcc_cell() {
  optimet::Matrix<optimet::t_real> result = optimet::Matrix<optimet::t_real>::Ones(3, 3) * 0.5;
  result.diagonal().fill(0);
  return result;
}
}

int main(int argc, char *const argv[]) {
  using namespace optimet;
#if defined(OPTIMET_MPI) && !defined(OPTIMET_JUST_DO_SERIAL)
  mpi::init(argc, const_cast<const char **>(argv));
#endif

  auto const iterations = 50;
  auto const nobjects = 100;
  auto const radius = 0.25;
  auto const nMax = 10;
  ElectroMagnetic const elmag{13.1, 1.0};
  auto const length = (radius + 0.5) * default_length();
  Scatterer const scatterer = {{0, 0, 0}, elmag, radius * default_length(), nMax};

  // setup geometry
  auto const cell = fcc_cell();
  t_int n = std::pow(nobjects, 1e0 / 3e0);
  auto range = std::make_tuple(n, n, n);
  if(n * n * n < nobjects)
    ++std::get<0>(range);
  if((n + 1) * n * n < nobjects)
    ++std::get<1>(range);
  auto geometry = std::make_shared<Geometry>();
  n = 0;
  for(int i(0); i < std::get<0>(range); ++i)
    for(int j(0); j < std::get<1>(range); ++j)
      for(int k(0); k < std::get<2>(range); ++k, ++n) {
        if(n == nobjects)
          break;
        Eigen::Matrix<t_real, 3, 1> pos = cell * Eigen::Matrix<t_real, 3, 1>(i, j, k) * length;
        geometry->pushObject({Tools::toSpherical({pos(0), pos(1), pos(2)}), scatterer.elmag,
                              scatterer.radius, scatterer.nMax});
      }

  // Create excitation
  Spherical<t_real> const vKinc{2 * consPi / default_wavelength(), 90 * consPi / 180.0,
                                90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto const excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, scatterer.nMax);
  excitation->populate();
  geometry->update(excitation);

#if defined(OPTIMET_MPI) && !defined(OPTIMET_JUST_DO_SERIAL)
  mpi::Communicator const world;
  auto const subdiagonals = std::max<int>(1, geometry->objects.size() / 2 - 2);
  mpi::FastMatrixMultiply const fmm(geometry->bground, excitation->wavenumber(), geometry->objects,
                                    subdiagonals, world);
#else
  FastMatrixMultiply const fmm(geometry->bground, excitation->wavenumber(), geometry->objects);
#endif

  Vector<t_complex> const input = Vector<t_complex>::Random(fmm.cols());

  Vector<t_complex> result(fmm.cols());
  auto start = std::chrono::high_resolution_clock::now();
  for(int i(0); i < iterations; ++i)
    fmm(input, result);
  auto end = std::chrono::high_resolution_clock::now();

  auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
#if defined(OPTIMET_MPI) && !defined(OPTIMET_JUST_DO_SERIAL)
  auto proc_max = world.all_reduce(elapsed_seconds.count(), MPI_MAX);
#else
  auto proc_max = elapsed_seconds.count();
#endif

#if defined(OPTIMET_MPI) && !defined(OPTIMET_JUST_DO_SERIAL)
  if(world.is_root())
#endif
    std::cout << "Timing: " << proc_max / iterations << "s\n";

#if defined(OPTIMET_MPI) && !defined(OPTIMET_JUST_DO_SERIAL)
  mpi::finalize();
#endif
  return 0;
}
