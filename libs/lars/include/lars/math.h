#pragma once

namespace lars {
  
  template<class T> inline constexpr T pow(const T base, const unsigned exponent) {
    return (exponent == 0) ? 1 : (exponent % 2 == 0) ? pow(base, exponent/2)*pow(base, exponent/2) : base * pow(base, (exponent-1)/2) * pow(base, (exponent-1)/2);
  }

  template <typename T> static int sign(T val) {
    return (T(0) < val) - (val < T(0));
  }
  

}