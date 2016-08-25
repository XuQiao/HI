#ifndef CORRELATIONS_CORRELATOR_HH
#define CORRELATIONS_CORRELATOR_HH
/**
 * @file   correlations/Correlator.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  Base class for correlators
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
#include <correlations/Result.hh>
#include <correlations/Types.hh>
#include <iostream>

namespace correlations {
  //____________________________________________________________________
  /**
   * Base class for calculating the correlations.
   *
   @code
   correlations::QVector q;
   correlations::RealVector phis;
   correlations::RealVector weights;
   correlations::Result  r;
   correlations::Correlator c = ...;
   correlations::HarmonicVector h(n);
   
   while (moreEvents) {
   
     while (moreObservations) {
       correlations::Real phi    = NextObservation();
       correlations::Real Weight = GetWeight(phi);
   
       phis.push_back(phi);
       weights.push_back(weight);
       q.fill(phi, weight);
     }
   
     r += c.calculate(h);
   }
   std::cout << r.eval() << std::endl;
   @endcode
   * @headerfile ""  <correlations/Correlator.hh>
   *
   */
  struct Correlator
  {
    enum {
      /** Value to mark ignored harmonics */
      kInvalidN = 0x7FFF
    };
    /**
     * Calculate up to the 8-particle correlator with the specified harmonics.
     *
     * @param h1 Harmonic 1
     * @param h2 Harmonic 2
     * @param h3 Harmonic 3
     * @param h4 Harmonic 4
     * @param h5 Harmonic 5
     * @param h6 Harmonic 6
     * @param h7 Harmonic 7
     * @param h8 Harmonic 8
     *
     * @return The correlation
     */
    Result calculate(const Harmonic h1,
		     const Harmonic h2=kInvalidN,
		     const Harmonic h3=kInvalidN,
		     const Harmonic h4=kInvalidN,
		     const Harmonic h5=kInvalidN,
		     const Harmonic h6=kInvalidN,
		     const Harmonic h7=kInvalidN,
		     const Harmonic h8=kInvalidN) const
    {
      const Harmonic hs[] = { h1, h2, h3, h4, h5, h6, h7, h8 };
      HarmonicVector hv;
      for (Size i = 0; i < 8; i++) {
	if (hs[i] == kInvalidN) break;
	hv.push_back(hs[i]);
      }
      return calculate(hv);
    }
    /**
     * Calculate the multi-particle correlation.
     *
     * @param h Vector of harmonics.  The correlator size is the size
     * of this vector.
     *
     * @return The correlation
     */
    Result calculate(const HarmonicVector& h) const
    {
      return calculate(h.size(), h);
    }
    /**
     * Calculate up the multi-particle  correlator.
     *
     * @param size How many particles to correlate
     * @param h    Harmonics of each term
     *
     * @return The correlation
     */
    Result calculate(const Size size, const HarmonicVector& h) const
    {
#if 0
      if (size > h.size()) {
	std::cerr << "Not enough harmonics (" << h.size() << ") to do "
		  << size << "-particle correlations" << std::endl;
	return Result();
      }
#endif
      if (size > _maxFixed) {
	// If we're above our fixed cut, or above 8, use the generic
	// algorithms
	return cN(size, h);
      }

      // Now, use specific implementations
      switch (size) {
      case 1: return c1(h);  break;
      case 2: return c2(h);  break;
      case 3: return c3(h);  break;
      case 4: return c4(h);  break;
      case 5: return c5(h);  break;
      case 6: return c6(h);  break;
      case 7: return c7(h);  break;
      case 8: return c8(h);  break;
      }
      std::cerr << "Number of correlators too big:" << size << std::endl;
      return Result();
    }
    /**
     * @return Name of the correlator
     */
    virtual const char* name() const { return "Correlator"; }
    /**
     * Destructor
     */
    virtual ~Correlator() {}
  protected:
    Correlator(Size maxFixed=8) : _maxFixed(maxFixed > 8 ? 8 : maxFixed) {}
    /**
     * Calculate the @a n particle correlation using harmonics @a h.
     *
     * @f[
     * C\{n\} = \langle\exp[i(\sum_j^n h_j\phi_j)]\rangle
     * @f]
     *
     * @param n How many particles to correlate
     * @param h @f$=h_1,\ldots,h_n@f$ Harmonic of each term
     *
     * @return The correlator and the summed weights
     */
    virtual Result cN(const Size n, const HarmonicVector& h) const = 0;
    /**
     * Calculate 1-particle correlator.
     *
     * @f[
     * C\{1\} = \langle\exp[i(\sum_j^1 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonic @f$ h_j@f$
     *
     * @return @f$ C\{1\}@f$
     */
    virtual Result c1(const HarmonicVector& h) const = 0;
    /**
     * Calculate 2-particle correlator.
     * @f[
     * C\{2\} = \langle\exp[i(\sum_j^2 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$ C\{2\}@f$
     */
    virtual Result c2(const HarmonicVector& h) const = 0;
    /**
     * Calculate 3-particle correlator.
     * @f[
     * C\{3\} = \langle\exp[i(\sum_j^3 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$QC\{3\}@f$
     */
    virtual Result c3(const HarmonicVector& h) const = 0;
    /**
     * Calculate 4-particle correlator.
     * @f[
     * C\{4\} = \langle\exp[i(\sum_j^4 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$QC\{4\}@f$
     */
    virtual Result c4(const HarmonicVector& h) const = 0;
    /**
     * Calculate 5-particle correlator.
     * @f[
     * C\{5\} = \langle\exp[i(\sum_j^5 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$QC\{5\}@f$
     */
    virtual Result c5(const HarmonicVector& h) const =0;
    /**
     * Calculate 6-particle correlator.
     * @f[
     * C\{6\} = \langle\exp[i(\sum_j^6 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$QC\{6\}@f$
     */
    virtual Result c6(const HarmonicVector& h) const =0;
    /**
     * Calculate 7-particle correlator.
     * @f[
     * C\{7\} = \langle\exp[i(\sum_j^7 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$QC\{7\}@f$
     */
    virtual Result c7(const HarmonicVector& h) const =0;
    /**
     * Calculate 8-particle correlator.
     * @f[
     * C\{8\} = \langle\exp[i(\sum_j^8 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$QC\{8\}@f$
     */
    virtual Result c8(const HarmonicVector& h) const =0;

    Size _maxFixed;
  };
}
#endif
// Local Variables:
//  mode: C++
// End:
