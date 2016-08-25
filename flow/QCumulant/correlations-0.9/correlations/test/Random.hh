#ifndef CORRELATIONS_TEST_RANDOM_H
#define CORRELATIONS_TEST_RANDOM_H
/**
 * @file   correlations/test/Random.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  Random number generator
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
#include <cstdlib>
#include <ctime>

namespace correlations {
  namespace test {
    /**
     * Random number generator.
     *
     * @headerfile "" <correlations/test/Random.hh>
     */
    struct Random
    {
      static void seed(unsigned int seed)
      {
	srand(seed == 0 ? time(0) : seed);
      }
      // using Real;
      // using Size;
      static Real asReal(Real min, Real max)
      {
        return min + (max-min) * Real(rand()) / RAND_MAX;
      }
      static Harmonic asHarmonic(Harmonic min, Harmonic max)
      {
    	return min + Harmonic((max-min) * Real(rand()) / RAND_MAX);
      }
      static Size asSize(Size min, Size max)
      {
	return min + Size((max-min) * Real(rand()) / RAND_MAX);
      }
    };
  }
}
#endif
// Local Variables:
//  mode: C++
// End:
