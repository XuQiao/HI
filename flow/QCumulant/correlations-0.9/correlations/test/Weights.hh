#ifndef WEIGHTS_H
#define WEIGHTS_H
/**
 * @file   correlations/test/Weights.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  Particle weights
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
#include <correlations/Types.hh>
#include <correlations/test/Random.hh>

namespace correlations {
  namespace test {
    /**
     * A structure that hold test weights. 
     *
     * The weights are random numbers between @f$ 1 \pm 0.144@@f$ and
     * depend on @f$ \phi@f$.  There are a total of 72 bins between
     * @f$ 0@f$ and @f$ 2\pi@f$.
     *
     * @headerfile "" <correlations/test/Weights.hh>
     */
    struct Weights
    {
      /**
       * Constructor
       *
       * @param n How many bins to make
       */
      Weights(Size n=72)
	: _value(n)
      {
	for (Size i = 0; i < n; i++)
	  _value[i] = 1 + Random::asReal(-0.144, 0.144);
      }
      /**
       * Find the weight corresponding to @a phi
       *
       * @param phi Azimuthal angle
       *
       * @return Weight corresponding to phi
       */
      Real eval(Real phi)
      {
	Size bin = phi / 2 / M_PI * _value.size();
	return _value[bin];
	// return 1;
      }
      /** Cache */
      RealVector _value;
    };
  }
}
#endif
// Local Variables:
//  mode: C++
// End:
