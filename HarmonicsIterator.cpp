#include "Types.h"
#include "HarmonicsIterator.h"

namespace optimet {
void HarmonicsIterator::throw_on_invalid_nm(t_uint n, t_int m) {
  if(not valid_nm(n, m))
    throw std::runtime_error("abs(m) > n");
}

HarmonicsIterator& HarmonicsIterator::principal(t_uint n) {
  principal_ = n;
  if(not valid_nm(n, m()))
    secondary_ = -static_cast<t_int>(n);
  return *this;
}

HarmonicsIterator& HarmonicsIterator::secondary(t_int m) {
  secondary_ = m;
  if(not valid_nm(n(), m))
    principal_ = static_cast<t_uint>(std::abs(m));
  return *this;
}

HarmonicsIterator HarmonicsIterator::operator++(int) {
  HarmonicsIterator const c(*this);
  operator++();
  return c;
}
HarmonicsIterator& HarmonicsIterator::operator++() {
  --secondary_;
  if(not valid_nm(principal_, secondary_)) {
    ++principal_;
    secondary_ = static_cast<t_int>(principal_);
  }
  return *this;
}
HarmonicsIterator HarmonicsIterator::operator--(int) {
  HarmonicsIterator const c(*this);
  operator--();
  return c;
}
HarmonicsIterator& HarmonicsIterator::operator--() {
  // Absolute limit when decrementing
  if(principal_ == std::numeric_limits<t_uint>::max())
    return *this;

  ++secondary_;
  if(not valid_nm(principal_, secondary_)) {
    --principal_;
    secondary_ = -static_cast<t_int>(principal_);
  }
  return *this;
}

} /* optimet */
