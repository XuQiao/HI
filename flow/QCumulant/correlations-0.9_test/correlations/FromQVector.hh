#ifndef CORRELATIONS_FROMQVECTOR_HH
#define CORRELATIONS_FROMQVECTOR_HH
/**
 * @file   correlations/FromQVector.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  Base class for cumulant code
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
#include <correlations/QVector.hh>
#include <correlations/Correlator.hh>


namespace correlations {
  //____________________________________________________________________
  /**
   * Base class to calculate Cumulants of a given Q vector.
   *
   @code
   correlations::Result  r;
   correlations::HarmonicVector h(n);
   correlations::QVector q(h);
   correlations::FromQVector c(q);
   
   while (moreEvents) {
     q.reset();
   
     while (moreObservations) {
       correlations::Real phi    = NextObservation();
       correlations::Real Weight = GetWeight(phi);
   
       q.fill(phi, weight);
     }
   
     r += c.calculate(h);
   }
   std::cout << r.eval() << std::endl;
   @endcode
   * @headerfile ""  <correlations/FromQVector.hh>
   */
  struct FromQVector : public Correlator
  {
    virtual ~FromQVector() {}
  protected:
    /**
     * Constructor
     *
     * @param q Q vector to use
     */
    FromQVector(QVector& q) : Correlator(0), _q(q) {}
    /**
     * @{
     * @name Drivers
     */
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
      static HarmonicVector null(0);
      if (null.size() <= 0) null.resize(_q.maxN());
      static bool filled = false;
      if (!filled) {
        std::fill(null.begin(), null.end(), 0);
        filled = true;
      }
      return Result(ucN(n, h), ucN(n, null).real());
    }
    /**
     * Calculate 1-particle correlation.
     * @f[
     * QC\{1\} = \langle\exp[i(h_1\phi_1)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$ QC\{1\}@f$
     */
    virtual Result c1(const HarmonicVector& h) const
    {
      return Result(uc1(h[0]), uc1(0).real());
    }
    /**
     * Generic two-particle correlation.
     * @f[
     * QC\{2\} = \langle\exp[i(\sum_j^2 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$ QC\{2\}@f$
     */
    virtual Result c2(const HarmonicVector& h) const
    {
      return Result(uc2(h[0], h[1]), uc2(0, 0).real());
    }
    /**
     * Generic 3-particle correlation.
     * @f[
     * QC\{3\} = \langle\exp[i(\sum_j^3 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$ QC\{3\}@f$
     */
    virtual Result c3(const HarmonicVector& h) const
    {
      return Result(uc3(h[0], h[1], h[2]), uc3(0,0,0).real());
    }
    /**
     * Generic 4-particle correlation.
     * @f[
     * QC\{4\} = \langle\exp[i(\sum_j^4 h_j\phi_j)]\rangle
     * @f]
     *
     * @param h Harmonics @f$ h_j@f$
     *
     * @return @f$ QC\{4\}@f$
     */
    virtual Result c4(const HarmonicVector& h) const
    {
      return Result(uc4(h[0],h[1],h[2],h[3]), uc4(0,0,0,0).real());
    }
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
    virtual Result c5(const HarmonicVector& h) const
    {
      return Result(uc5(h[0], h[1], h[2], h[3], h[4]),
                    uc5(0, 0, 0, 0, 0).real());
    }
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
    virtual Result c6(const HarmonicVector& h) const
    {
      return Result(uc6(h[0],h[1],h[2],h[3],h[4],h[5]),
                    uc6(0,   0,   0,   0,   0,   0).real());
    }
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
    virtual Result c7(const HarmonicVector& h) const
    {
      std::cout << "Calculating C7" << std::endl;
      return Result(uc7(h[0],h[1],h[2],h[3],h[4],h[5],h[6]),
                    uc7(0,   0,   0,   0,   0,   0,   0).real());
    }
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
    virtual Result c8(const HarmonicVector& h) const
    {
      std::cout << "Calculating C8" << std::endl;
      return Result(uc8(h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7]),
                    uc8(0,   0,   0,   0,   0,   0,   0,   0).real());
    }
    /* @} */
    /**
     * @{
     * @name Abstract interface for fixed calculations
     */
    /**
     * Calculate the @f$ n@f$$-paricle correlator using harmonics
     * @f$ h = (h_1,\ldots,h_2)@f$
     *
     * @param n The number of particles to correlate
     * @param h Vector of harmonics.  @f$ n@f$ is set to the si
     *
     * @return @f$ \langle\exp[i(\sum_j^n h_j\phi_j)]\rangle@f$
     */
    virtual Complex ucN(const Size n, const HarmonicVector& h) const = 0;
    /**
     * Generic 1-particle correlation.
     * @f[
     * QC\{1\} = \langle\exp[i(\sum_j^1 h_j\phi_j)]\rangle
     * @f]
     *
     * @param n1 Harmonics @f$ h_1@f$
     *
     * @return @f$ QC\{1\}@f$
     */
    virtual Complex uc1(const Harmonic n1) const
    {
      return _q(n1, 1);
    }
    /**
     * Generic 2-particle correlation.
     * @f[
     * QC\{2\} = \langle\exp[i(\sum_j^2 h_j\phi_j)]\rangle
     * @f]
     *
     * @param n1 1st Harmonic
     * @param n2 2nd Harmonic
     *
     * @return the correlator
     */
    virtual Complex uc2(const Harmonic n1, const Harmonic n2) const
    {
      HarmonicVector h(2);
      h[0] = n1; h[1] = n2;
      return ucN(h.size(), h);
    }
    /**
     * Generic 3-particle correlation.
     * @f[
     * QC\{3\} = \langle\exp[i(\sum_j^3 h_j\phi_j)]\rangle
     * @f]
     *
     * @param n1 1st Harmonic
     * @param n2 2nd Harmonic
     * @param n3 3rd Harmonic
     *
     * @return the correlator
     */
    virtual Complex uc3(const Harmonic n1,
                       const Harmonic n2,
                       const Harmonic n3) const
    {
      HarmonicVector h(3);
      h[0] = n1; h[1] = n2; h[2] = n3;
      return ucN(h.size(), h);
    }
    /**
     * Generic 4-particle correlation.
     * @f[
     * QC\{4\} = \langle\exp[i(\sum_j^4 h_j\phi_j)]\rangle
     * @f]
     *
     *
     * @param n1 1st Harmonic
     * @param n2 2nd Harmonic
     * @param n3 3rd Harmonic
     * @param n4 4th Harmonic
     *
     * @return the correlator
     */
    virtual Complex uc4(const Harmonic n1,
                        const Harmonic n2,
                       const Harmonic n3,
                       const Harmonic n4) const
    {
      HarmonicVector h(4);
      h[0] = n1; h[1] = n2; h[2] = n3; h[3] = n4;
      return ucN(h.size(), h);
    }
    /**
     * Generic 5-particle correlation.
     * @f[
     * QC\{5\} = \langle\exp[i(\sum_j^5 h_j\phi_j)]\rangle
     * @f]
     *
     * @param n1 1st Harmonic
     * @param n2 2nd Harmonic
     * @param n3 3rd Harmonic
     * @param n4 4th Harmonic
     * @param n5 5th Harmonic
     *
     * @return The correlator
     */
    virtual Complex uc5(const Harmonic n1,
                       const Harmonic n2,
                       const Harmonic n3,
                       const Harmonic n4,
                       const Harmonic n5) const
    {
      HarmonicVector h(5);
      h[0] = n1; h[1] = n2; h[2] = n3; h[3] = n4;
      h[4] = n5;
      return ucN(h.size(), h);
    }
    /**
     * Calculate 6-particle correlator
     * @f[
     * C\{6\} = \langle\exp[i(\sum_j^6 h_j\phi_j)]\rangle
     * @f]
     *
     * @param n1 1st Harmonic
     * @param n2 2nd Harmonic
     * @param n3 3rd Harmonic
     * @param n4 4th Harmonic
     * @param n5 5th Harmonic
     * @param n6 6th Harmonic
     *
     * @return @f$QC\{6\}@f$
     */
    virtual Complex uc6(const Harmonic n1,
                       const Harmonic n2,
                       const Harmonic n3,
                       const Harmonic n4,
                       const Harmonic n5,
                       const Harmonic n6) const
    {
      HarmonicVector h(6);
      h[0] = n1; h[1] = n2; h[2] = n3; h[3] = n4;
      h[4] = n5; h[5] = n6; // h[6] = n7; h[7] = n8;
      return ucN(h.size(), h);
    }
    /**
     * Calculate 7-particle correlator
     * @f[
     * C\{7\} = \langle\exp[i(\sum_j^7 h_j\phi_j)]\rangle
     * @f]
     *
     * @param n1 1st Harmonic
     * @param n2 2nd Harmonic
     * @param n3 3rd Harmonic
     * @param n4 4th Harmonic
     * @param n5 5th Harmonic
     * @param n6 6th Harmonic
     * @param n7 7th Harmonic
     *
     * @return @f$QC\{7\}@f$
     */
    virtual Complex uc7(const Harmonic n1,
                       const Harmonic n2,
                       const Harmonic n3,
                       const Harmonic n4,
                       const Harmonic n5,
                       const Harmonic n6,
                       const Harmonic n7) const
    {
      HarmonicVector h(7);
      h[0] = n1; h[1] = n2; h[2] = n3; h[3] = n4;
      h[4] = n5; h[5] = n6; h[6] = n7; // h[7] = n8;
      return ucN(h.size(), h);
    }
    /**
     * Calculate 8-particle correlator
     * @f[
     * C\{8\} = \langle\exp[i(\sum_j^8 h_j\phi_j)]\rangle
     * @f]
     *
     * @param n1 1st Harmonic
     * @param n2 2nd Harmonic
     * @param n3 3rd Harmonic
     * @param n4 4th Harmonic
     * @param n5 5th Harmonic
     * @param n6 6th Harmonic
     * @param n7 7th Harmonic
     * @param n8 8th Harmonic
     *
     * @return @f$QC\{8\}@f$
     */
    virtual Complex uc8(const Harmonic n1,
                       const Harmonic n2,
                       const Harmonic n3,
                       const Harmonic n4,
                       const Harmonic n5,
                       const Harmonic n6,
                       const Harmonic n7,
                       const Harmonic n8) const
    {
      HarmonicVector h(8);
      h[0] = n1; h[1] = n2; h[2] = n3; h[3] = n4;
      h[4] = n5; h[5] = n6; h[6] = n7; h[7] = n8;
      return ucN(h.size(), h);
    }
    /* @} */
    QVector& _q;
  };
}
#endif
// Local Variables:
//  mode: C++
// End:
