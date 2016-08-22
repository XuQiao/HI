#ifndef CORRELATIONS_TYPES_HH
#define CORRELATIONS_TYPES_HH
/**
 * @file   correlations/Types.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  Declarations of types
 */
/*
 * Multi-particle correlations 
 * Copyright (C) 2013 K.Gulbrandsen, A.Bilandzic, C.H. Christensen.
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses.
 */
#include <vector>
#include <complex>

/**
 * Namespace for correlations code
 */
namespace correlations {
  /**
   * Type of Real.
   *
   * Redefining to this to be long double say would make all of the code
   * use that precsion.
   */
  typedef double Real;
  /** Type of real array */
  typedef std::vector<Real> RealVector;
  /** Type defintion of the kind of complex we use */
  typedef std::complex<Real> Complex;
  /** Type definition of our internal storage */
  typedef std::vector<Complex> ComplexVector;
  /** Type of harmonics */
  typedef short Harmonic;
  /**
   * Type definition of vector of harmonics
   *
   * @code
   * correlations::Result  r;
   * correlations::HarmonicVector h(n);
   * for (size_t i = 0; i < n; i++) h[i] = ...;
   * correlations::Correlator c = ...;n
   *
   * while (moreEvents) {
   *
   *   while (moreObservations) {
   *     ...
   *   }
   *
   *   r += c.calculate(h);
   * }
   * std::cout << r.eval() << std::endl;
   * @endcode
   */
  typedef std::vector<Harmonic> HarmonicVector;
  /** Type of Sizes */
  typedef unsigned short Size;
  /** Type of Powers */
  typedef unsigned short Power;

}
#endif
// Local Variables:
//   mode: C++
// End:
