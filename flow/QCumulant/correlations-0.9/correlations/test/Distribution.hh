#ifndef CORRELATIONS_TEST_DISTRIBUTION_H
#define CORRELATIONS_TEST_DISTRIBUTION_H
/**
 * @file   correlations/test/Distribution.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  A phi distribution
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
#include <correlations/test/Random.hh>
#include <correlations/Types.hh>
#include <cmath>
#include <vector>
#include <algorithm>

namespace correlations {
  namespace test {
    //====================================================================
    /**
     * The @f$\phi@f$ distribution.
     *
     * @f[
     *  f(\phi) = 1 + 2\sum_{i=1}^{6} p_i\cos(i(\phi-\Phi_R))
     * @f]
     *
     * The indefinite integral is given by
     * @f[
     * \int d\phi f(\phi) = 2\pi + 
     *    2 \sum_{i=1}^{6} p_i \frac{\sin(i(\phi-\Phi_R))}{i}
     * @f]
     *
     * or for the definte integral over @f$[a,b]@f$
     * @f[
     * \int_a^b f(\phi) = (b-a) + 2 \sum_{i=1}^6
     *   p_i \frac{\sin(i(b-\Phi_R))-\sin(i(a-\Phi_R))}{i}
     * @f]
     *
     * @headerfile "" <correlations/test/Distribution.hh>
     */
    struct Distribution
    {
      static Real func(Real* xx,
		       Real* pp)
      {
	Real phi = xx[0];
	Real psi = pp[0];
	Real ret = 1;
	for (Size i = 1; i < 7; i++)
	  ret += 2 * pp[i] * cos(i * (phi-psi));
	return ret;
      }

      /**
       * Constructor
       */
      Distribution()
	: _npx(100),
	  _integral(_npx+1),
	  _alpha(_npx),
	  _beta(_npx),
	  _gamma(_npx)
      {
	_v[0] = 0;    // psi
	_v[1] = 0.01; // v_1
	_v[2] = 0.05; // v_2
	_v[3] = 0.03; // v_3
	_v[4] = 0.02; // v_4
	_v[5] = 0.01; // v_5
	_v[6] = 0.01; // v_6;
      }
      /**
       * Set-up for sampling
       */
      void setup(Real psi)
      {
	_v[0] = psi;
	Real dx    = 2*M_PI / _npx;
	_integral[0] = 0;
	// Calculate cumulated integral in bins
	for (Size i = 0; i < _npx; i++) {
	  Real intg = integral(i*dx, (i+1)*dx);
	  _integral[i+1] = _integral[i] + intg;
	}
	// Normalize intergral
	Real full = _integral[_npx];
	for (Size i = 0; i <= _npx; i++)
	  _integral[i] /= full;
	// Approximate intergal r in each bin with a parabolla
	//
	//  x = alpha + beta * r + gamma * r * r;
	for (Size i = 0; i <= _npx; i++) {
	  Real r2 = _integral[i+1] - _integral[i];     // End point
	  Real r1 = integral(i*dx, (i+.5)*dx) / full;  // Mid point
	  Real r3 = 2 * r2 - 3 * r1;
	  _gamma[i] = 2 * (fabs(r3) > 1e-8) ? r3 / (dx*dx) : 0;
	  _beta[i]  = r2 / dx - _gamma[i] * dx / 2;
	  _alpha[i] = i * dx;
	}
      }
      /**
       * Evaluate the function
       *
       * @param phi Where to evaluate
       *
       * @return The function evaluated at @a phi
       */
      Real eval(Real phi) const
      {
	Real* pp = const_cast<Real*>(&(_v[0]));
	return func(&phi, pp);
      }
      /**
       * Evaluate the integral of the function from @a a to @a b
       *
       * @param a Lower limit
       * @param b Upper limit
       *
       * @return Integral from @a a to @a b
       */
      Real integral(Real a,
		    Real b) const
      {
	Real ret = (b-a);
	Real psi = _v[0];
	for (Size i = 1; i < 7; i++)
	  ret += 2 * _v[i] * (sin(i*(b-psi))-sin(i*(a-psi)));
	return ret;
      }
      /**
       * Draw a random value from the function considered a probability
       * distribution.
       *
       * @return A random number distributed according to the function
       */
      Real random() const
      {
	Real        rnd   = Random::asReal(0,1);
	const Real* first = &(_integral[0]);
	const Real* last  = first + _integral.size();
	const Real* ptr   = std::lower_bound(first, last, rnd);
	Size        bin   = 0;
	if (ptr != (last) && (*ptr == rnd)) bin = (ptr-first);
	else                                bin = (ptr-first-1);
	Real  rr    = rnd - _integral[bin];

	Real yy = 0;
	Real g  = _gamma[bin];
	Real b  = _beta[bin];
	Real a  = _alpha[bin];
	if (g != 0)
	  yy = (-b + sqrt(b*b+2*g*rr))/g;
	else
	  yy = rr / b;
	Real x = a + yy;
	return x;
      }
      /** Reaction plane and Flow strength */
      Real _v[7];
      /** Number of points to evaluate the function at */
      Size _npx;
      /** Integral of function in steps */
      RealVector _integral;
      /** Parameterization of integral */
      RealVector _alpha;
      /** Parameterization of integral */
      RealVector _beta;
      /** Parameterization of integral */
      RealVector _gamma;
    };
  }
}
#endif
// Local Variables:
//  mode: C++
// End:
