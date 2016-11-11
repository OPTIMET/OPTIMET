#ifndef OPTIMET_COEFFICIENTS_H
#define OPTIMET_COEFFICIENTS_H

#include "Types.h"
#include <cmath>

namespace optimet {
namespace coefficient {
template <class TYPE = t_real> TYPE a(t_uint n, t_int m) {
  t_uint const absm(std::abs(m));
  if(n < absm)
    return static_cast<TYPE>(0);
  return std::sqrt(static_cast<TYPE>((n + 1 + absm) * (n + 1 - absm)) /
                   static_cast<TYPE>((2 * n + 1) * (2 * n + 3)));
}

template <class TYPE = t_real> TYPE b(t_uint n, t_int m) {
  if(static_cast<t_uint>(std::abs(m)) > n)
    return static_cast<TYPE>(0);
  return (m >= 0 ? 1 : -1) * std::sqrt(static_cast<TYPE>((n - m - 1) * (n - m)) /
                                       static_cast<TYPE>((2 * n - 1) * (2 * n + 1)));
}

template <class TYPE = t_real> TYPE c(t_uint n, t_int m) {
  if(std::abs(m) > n)
    return 0;
  return (m < 0 ? -1 : 1) *
         std::sqrt((static_cast<t_int>(n) - m) * (static_cast<t_int>(n) + m + 1));
}
}
}
#endif
