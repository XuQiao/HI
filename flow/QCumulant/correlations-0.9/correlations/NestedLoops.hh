#ifndef CORRELATIONS_NESTEDLOOPS_HH
#define CORRELATIONS_NESTEDLOOPS_HH
/**
 * @file   correlations/NestedLoops.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  Concrete and base class for nested loops
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
#include <cmath>
#include <correlations/Correlator.hh>

namespace correlations {
  //____________________________________________________________________
  /**
   * Structure to calculate the correlator using nested loops
   *
   @code
   correlations::Result         r;
   correlations::RealVector     phis;
   correlations::RealVector     weights;
   correlations::NestedLoops    c(phis,weights);
   correlations::HarmonicVector h(n);
   
   while (moreEvents) {
   
     while (moreObservations) {
       correlations::Real phi    = NextObservation();
       correlations::Real weight = GetWeight(phi);
   
       phis.push_back(phi);
       weights.push_back(weight);
     }
   
     r += c.calculate(h);
   }
   std::cout << r.eval() << std::endl;
   @endcode
   * @headerfile ""  <correlations/NestedLoops.hh>
   */
  struct NestedLoops : public Correlator
  {
    /**
     * Constructor
     *
     * @param phis       Reference to phi array that will  be filled
     * @param weights    Reference to weight array that will be filled
     * @param useWeights Whether to use weights or not
     */
    NestedLoops(RealVector& phis,
		RealVector& weights,
		bool        useWeights=true)
      : Correlator(6),
	_useWeights(useWeights),
	_phis(phis),
	_weights(weights)
    {}
    /**
     * @return Name of the correlator
     */
    virtual const char* name() const { return "Explicit nested loops"; }
  protected:
    /**
     * Update the elements in @a p and @a w for order @a iorder and
     * particle @a z[iorder].  But first, the array @a z is checked to
     * see if @a z[iorder] is present - if so, we return false to
     * avoid auto-correlations
     *
     * @param iorder Current index into @a p, @a w, and @a z
     * @param z Indices of particles
     * @param p Particle @f$ \phi@f$ weighted by the harmonic
     * @param w Weight of particles
     * @param h Harmonic vector
     *
     * @return true if we should go on with this
     */
    bool store(Size iorder, Size* z, Real* p, Real* w,
	       const HarmonicVector& h) const
    {
      Size*  first = z;
      Size*  last  = z+iorder;
      Size   ipart = z[iorder];
      // Check for auto-correlation
      Size*  pp = first;
      while (pp < last) {
 	if (*pp == ipart) {
	  return false;
	}
	pp++;
      }
#if 0
      std::cout << std::endl;
      if (std::find(first, last, ipart) != last) {
	std::ostream_iterator<Size> o(std::cerr, " ");
	std::cerr << iorder << ": " << ipart << " found in ";
	std::copy(first, last, o);
	std::cerr << " - auto correlation" << std::endl;
	return false;
      }
#endif
      p[iorder] = h[iorder] * _phis[ipart];
      w[iorder] = _weights[ipart];
      return true;
    }

    /**
     * Calculate a single term
     *
     * @param l      Number of arguments in @a w and @a p
     * @param w      Array of weights
     * @param p      Array of phis multiplied by harmonic order
     * @param c      Complex number to add to
     * @param sumw   Sum of weights are calculated here
     *
     * @return The complex term to add
     */
    void term(const Size  l,
	      const Real* w,
	      const Real* p,
	      Complex&    c,
	      Real&       sumw) const
    {
      Real arg = 0;
      Real ww  = 1;
      for (Size i = 0; i < l; i++) {
	arg += p[i];
	if (_useWeights) ww *= w[i];
      }
      sumw += ww;
      c    += Complex(ww*cos(arg), ww*sin(arg));
    }
    /**
     * Calculate the @a n particle correlation using harmonics @a h
     *
     * @param n How many particles to correlate
     * @param h Harmonic of each term
     *
     * @return The correlator and the summed weights
     */
    virtual Result cN(const Size n, const HarmonicVector& h) const
    {
      switch (n) {
      case 1: return c1(h);   break;
      case 2: return c2(h);   break;
      case 3: return c3(h); break;
      case 4: return c4(h);  break;
      case 5: return c5(h);  break;
      case 6: return c6(h);  break;
      case 7: return c7(h);  break;
      case 8: return c8(h);  break;
      }
      std::cerr << "Number of correlators too big:" << n << std::endl;
      return Result();
    }
    /**
     * @{
     * @name Fixed implementations
     */
    /**
     * Calculate correlator for 1 particle
     *
     * @param h Harmonics
     *
     * @return Correlator
     */
    Result c1(const HarmonicVector& h) const
    {
      Complex r;
      Real sumw = 0;
      for (Size i = 0; i < _phis.size(); i++) {
	Real w = _weights[i];
	Real p = h[0] * _phis[i];
	term(1, &w, &p, r, sumw);
      }
      return Result(r, sumw);
    }
    enum {
      kInvalidIdx = 0x7FFF
    };
    /**
     * Calculate correlator for 2 particles.
     *
     * @param h Harmonics
     *
     * @return Correlator
     */
    Result c2(const HarmonicVector& h) const
    {
      Complex      r;
      Real         sumw = 0;
      Real         p[2];
      Real         w[2];
      Size         z[2] = { kInvalidIdx, kInvalidIdx };
      for (z[0] = 0; z[0] < _phis.size(); z[0]++) {
	if (!store(0,z,p,w,h)) continue;
	for (z[1] = 0; z[1] < _phis.size(); z[1]++) {
	  if (!store(1,z,p,w,h)) continue;
	  term(2, w, p, r, sumw);
	}
      }
      return Result(r, sumw);
    }
    /**
     * Calculate correlator for 3 particles.
     *
     * @param h Harmonics
     *
     * @return Correlator
     */
    Result c3(const HarmonicVector& h) const
    {
      Complex r;
      Real sumw = 0;
      Real p[3];
      Real w[3];
      Size z[3] = { kInvalidIdx, kInvalidIdx, kInvalidIdx };
      for (z[0] = 0; z[0] < _phis.size(); z[0]++) {
	if (!store(0,z,p,w,h)) continue;
	for (z[1] = 0; z[1] < _phis.size(); z[1]++) {
	  if (!store(1,z,p,w,h)) continue;
	  for (z[2] = 0; z[2] < _phis.size(); z[2]++) {
	    if (!store(2,z,p,w,h)) continue;
	    term(3, w, p, r, sumw);
	  }
	}
      }
      return Result(r, sumw);
    }
    /**
     * Calculate correlator for 4 particles.
     *
     * @param h Harmonics
     *
     * @return Correlator
     */
    Result c4(const HarmonicVector& h) const
    {
      Complex r;
      Real sumw = 0;
      Real p[4];
      Real w[4];
      Size z[4] = { kInvalidIdx, kInvalidIdx, kInvalidIdx, kInvalidIdx };
      for (z[0] = 0; z[0] < _phis.size(); z[0]++) {
	if (!store(0,z,p,w,h)) continue;
	for (z[1] = 0; z[1] < _phis.size(); z[1]++) {
	  if (!store(1,z,p,w,h)) continue;
	  for (z[2] = 0; z[2] < _phis.size(); z[2]++) {
	    if (!store(2,z,p,w,h)) continue;
	    for (z[3] = 0; z[3] < _phis.size(); z[3]++) {
	      if (!store(3,z,p,w,h)) continue;
	      term(4, w, p, r, sumw);
	    }
	  }
	}
      }
      return Result(r, sumw);
    }
    /**
     * Calculate 5-particle correlator.
     * @f[
     * C\{5\} = \langle\exp[i(\sum_z[1]^5 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$QC\{5\}@f$
     */
    virtual Result c5(const HarmonicVector& h) const
    {
      Complex r;
      Real sumw = 0;
      Real p[5];
      Real w[5];
      Size z[5] = { kInvalidIdx, kInvalidIdx, kInvalidIdx, kInvalidIdx,
		    kInvalidIdx };
      for (z[0] = 0; z[0] < _phis.size(); z[0]++) {
	store(0,z,p,w,h);
	for (z[1] = 0; z[1] < _phis.size(); z[1]++) {
	  if (!store(1,z,p,w,h)) continue;
	  for (z[2] = 0; z[2] < _phis.size(); z[2]++) {
	    if (!store(2,z,p,w,h)) continue;
	    for (z[3] = 0; z[3] < _phis.size(); z[3]++) {
	      if (!store(3,z,p,w,h)) continue;
	      for (z[4] = 0; z[4] < _phis.size(); z[4]++) {
		if (!store(4,z,p,w,h)) continue;
		term(5, w, p, r, sumw);
	      } // for z[4]
	    } // for z[3]
	  } // for z[2]
	} // for z[1]
      } // for z[0]
      return Result(r, sumw);
    }
    /**
     * Calculate 6-particle correlator.
     * @f[
     * C\{6\} = \langle\exp[i(\sum_z[1]^6 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$QC\{6\}@f$
     */
    virtual Result c6(const HarmonicVector& h) const
    {
      Complex r;
      Real sumw = 0;
      Real p[6];
      Real w[6];
      Size z[6] = { kInvalidIdx, kInvalidIdx, kInvalidIdx, kInvalidIdx,
		    kInvalidIdx, kInvalidIdx };
      for (z[0] = 0; z[0] < _phis.size(); z[0]++) {
	store(0,z,p,w,h);
	for (z[1] = 0; z[1] < _phis.size(); z[1]++) {
	  if (!store(1,z,p,w,h)) continue;
	  for (z[2] = 0; z[2] < _phis.size(); z[2]++) {
	    if (!store(2,z,p,w,h)) continue;
	    for (z[3] = 0; z[3] < _phis.size(); z[3]++) {
	      if (!store(3,z,p,w,h)) continue;
	      for (z[4] = 0; z[4] < _phis.size(); z[4]++) {
		if (!store(4,z,p,w,h)) continue;
		for (z[5] = 0; z[5] < _phis.size(); z[5]++) {
		  if (!store(5,z,p,w,h)) continue;
		  term(6, w, p, r, sumw);
		} // for z[5]
	      } // for z[4]
	    } // for z[3]
	  } // for z[2]
	} // for z[1]
      } // for z[0]
      return Result(r, sumw);
    }
    /**
     * Calculate 7-particle correlator.
     * @f[
     * C\{7\} = \langle\exp[i(\sum_z[1]^7 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$QC\{7\}@f$
     */
    virtual Result c7(const HarmonicVector& h) const
    {
      Complex r;
      Real sumw = 0;
      Real p[7];
      Real w[7];
      Size z[7] = { kInvalidIdx, kInvalidIdx, kInvalidIdx, kInvalidIdx,
		    kInvalidIdx, kInvalidIdx, kInvalidIdx };
      for (z[0] = 0; z[0] < _phis.size(); z[0]++) {
	store(0,z,p,w,h);
	for (z[1] = 0; z[1] < _phis.size(); z[1]++) {
	  if (!store(1,z,p,w,h)) continue;
	  for (z[2] = 0; z[2] < _phis.size(); z[2]++) {
	    if (!store(2,z,p,w,h)) continue;
	    for (z[3] = 0; z[3] < _phis.size(); z[3]++) {
	      if (!store(3,z,p,w,h)) continue;
	      for (z[4] = 0; z[4] < _phis.size(); z[4]++) {
		if (!store(4,z,p,w,h)) continue;
		for (z[5] = 0; z[5] < _phis.size(); z[5]++) {
		  if (!store(5,z,p,w,h)) continue;
		  for (z[6] = 0; z[6] < _phis.size(); z[6]++) {
		    if (!store(6,z,p,w,h)) continue;
		    term(7, w, p, r, sumw);
		  } // z[6]
		} // for z[5]
	      } // for z[4]
	    } // for z[3]
	  } // for z[2]
	} // for z[1]
      } // for z[0]
      return Result(r, sumw);
    }
    /**
     * Calculate 8-particle correlator.
     * @f[
     * C\{8\} = \langle\exp[i(\sum_z[1]^8 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$QC\{8\}@f$
     */
    virtual Result c8(const HarmonicVector& h) const
    {
      Complex r;
      Real sumw = 0;
      Real p[8];
      Real w[8];
      Size z[8] = { kInvalidIdx, kInvalidIdx, kInvalidIdx, kInvalidIdx,
		    kInvalidIdx, kInvalidIdx, kInvalidIdx, kInvalidIdx };
      for (z[0] = 0; z[0] < _phis.size(); z[0]++) {
	store(0,z,p,w,h);
	for (z[1] = 0; z[1] < _phis.size(); z[1]++) {
	  if (!store(1,z,p,w,h)) continue;
	  for (z[2] = 0; z[2] < _phis.size(); z[2]++) {
	    if (!store(2,z,p,w,h)) continue;
	    for (z[3] = 0; z[3] < _phis.size(); z[3]++) {
	      if (!store(3,z,p,w,h)) continue;
	      for (z[4] = 0; z[4] < _phis.size(); z[4]++) {
		if (!store(4,z,p,w,h)) continue;
		for (z[5] = 0; z[5] < _phis.size(); z[5]++) {
		  if (!store(5,z,p,w,h)) continue;
		  for (z[6] = 0; z[6] < _phis.size(); z[6]++) {
		    if (!store(6,z,p,w,h)) continue;
		    for (z[7] = 0; z[7] < _phis.size(); z[7]++) {
		      if (!store(7,z,p,w,h)) continue;
		      term(8, w, p, r, sumw);
		    } // for z[7]
		  } // for z[6]
		} // for z[5]
	      } // for z[4]
	    } // for z[3]
	  } // for z[2]
	} // for z[1]
      } // for z[0]
      return Result(r, sumw);
    }
    /** @} */
    bool        _useWeights;
    RealVector& _phis;
    RealVector& _weights;
  };
}
#endif
// Local Variables:
//  mode: C++
// End:
