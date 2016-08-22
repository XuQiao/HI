#ifndef CORRELATIONS_QVECTOR_HH
#define CORRELATIONS_QVECTOR_HH
/**
 * @file   correlations/QVector.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  Container of the Q-vector
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
#include <cmath>
#include <iostream>

namespace correlations {
  //____________________________________________________________________
  /**
   * A structure to hold the @f$Q@f$ vector.
   *
   * It provides member functions for reseting, filling and accessing
   * the individual @f$Q@f$ vector components.
   *
   @code
   correlations::QVector q(3,2);
   
   // Fill vector
   q.reset();
   for (Int_t i = 0; i < nPhi; i++) q.fill(phi[i], weight[i])
   
   // Show components
   for (short n=-2, n <= 2; n++)
     for (unsigned short p = 0; p < 2; p++)
       std::cout << "Q{" << n << "," << p << "}=" << q(n,p) << std::endl;
   @endcode
   * @headerfile ""  <correlations/QVector.hh>
   */
  struct QVector
  {
    /**
     * Constructor - creates a QVector with the specified maximum
     * harmonic order and maximum power of the harmonics
     *
     * @param mN         Maximum harmonic order
     * @param mP         Maximum power
     * @param useWeights Whether to use weights or not
     *
     */
    QVector():_maxN(0), _maxP(0), _useWeights(false), _q(0){}
    QVector(Size mN, Size mP, bool useWeights)
      : _maxN(mN), _maxP(mP), _useWeights(useWeights), _q(0)
    {
      resize(mN, mP);
    }
    /** 
     * Constructor 
     * 
     * @param h Harmonics vector
     * @param useWeights Whether to use weights
     */
    QVector(const HarmonicVector& h, bool useWeights) 
      : _maxN(0), _maxP(0), _useWeights(useWeights), _q(0)
    {
      resize(h);
    }
    /**
     * Resize the QVector - note this clears all content
     *
     * @param mN New maximum harmonic
     * @param mP New maximum power
     */
    void resize(Size mN, Size mP)
    {
      _maxN = mN;
      _maxP = mP;
      // We store both the positive and the negative harmonics for
      // efficiency reasons. So we need to index from -fMaxN-1 to fMaxN,
      // which gives us a total of 2*fMaxN+1 places for N.
      _q.resize((2*_maxN+1)*(_maxP+1));
    }
    /**
     * Resize the Q-vector to accommodate the harmonics specified in
     * the input harmonic vector. The max power is set to the size of
     * the input vector, and maximum harmonic is set the the sum of
     * the absolute value of each element in the vector.
     *
     * @param h Harmonic vector to resize to
     */
    void resize(const HarmonicVector& h)
    {
      _maxP = h.size();
      _maxN = 0;
      for (HarmonicVector::const_iterator i = h.begin(); i != h.end(); ++i)
	// Older GCC does not have overload for integer types
	_maxN += short(std::abs(float(*i)));
      resize(_maxN, _maxP);
    }
    /**
     * Reset all @f$Q@f$ vector component to @f$0+i0@f$
     *
     */
    void reset(const Complex& q=Complex(0,0))
    {
      std::fill(_q.begin(), _q.end(), q);
    }
    /**
     * @brief Get @f$Q_{n,p}@f$
     *
     * Get the @f$Q_{n,p} = \sum_{i=1}^{M} w_i^pe^{in\phi_i}@f$.
     * Note that @f$Q_{-n,p} = Q_{n,p}^*@f$
     *
     * @param n The harmonic
     * @param p The power
     *
     * @return @f$ Q_{n,p}@f$ or @f$ Q^*_{-n,p}@f$ if @f$n<0@f$
     */
    const Complex& operator()(Harmonic n, Power p) const
    {
      return _q[index(n,p)];
    }
    /**
    * @brief Get @f$Q_{n,p}@f$
    *
    * Get the @f$Q_{n,p} = \sum_{i=1}^{M} w_i^pe^{in\phi_i}@f$.
    * Note that @f$Q_{-n,p} = Q_{n,p}^*@f$
    *
    * @param n The harmonic
    * @param p The power
    *
    * @return @f$ Q_{n,p}@f$ or @f$ Q^*_{-n,p}@f$ if @f$n<0@f$
    */
    Complex& operator()(Harmonic n, Power p)
    {
      return _q[index(n,p)];
    }
    /**
     * Fill in an obersvation
     *
     * @param phi     Phi observation.
     * @param weight  Weight of this observation.
     */
    void fill(Real phi, Real weight)
    {
      std::vector<Real> cs(_maxN+1);
      std::vector<Real> sc(_maxN+1);
      for (Harmonic n = 0; n <= _maxN; n++) {
	cs[n] = cos(n * phi);
	sc[n] = sin(n * phi);
      }
      for (Power p=0; p <= _maxP; p++) {
	Real w = (!_useWeights || p == 0 ? 1 : pow(weight,p));
	for (Harmonic n = 0; n <= _maxN; n++) {
	  this->operator()(n,p) += Complex(w * cs[n], w * sc[n]);
	  if (n == 0) continue;
	  // Update the conjugate as well
	  this->operator()(-n,p) += Complex(w * cs[n], -w * sc[n]);
	}
      }
    }
    /**
     * Get the maximum harmonic (minus 1)
     *
     * @return Maximum harmonic minus 1
     */
    Size maxN() const { return _maxN; }
    /**
     * Get the maximum power (minus 1)
     *
     * @return Maximum power (minus 1)
     */
    Power maxP() const { return _maxP; }

    bool verify() const
    {
      bool ret = true;
      for (Harmonic i = -_maxN; i <= _maxN; i++) {
	if (i == 0) continue;
	for (Power j = 0; j <= _maxP; j++) {
	  const Complex& c1 = (*this)(i,j);
	  const Complex& c2 = (*this)(-i,j);
	  if (c1 != std::conj(c2)) {
	    std::cerr << "(" << i << "," << j << ") [" << index(i,j)
		      << "]=" << c1 << " is no the conjugate of ("
		      << -i << "," << j << ") [" << index(-i,j)
		      << "]= " << c2 << std::endl;
	    ret = false;
	  }
	}
      }
      return ret;
    }

    void print() const {
      if (!verify()) return;
      for (Harmonic i = -_maxN; i <= _maxN; i++) {
	std::cout << "Harmonic # " << i << std::endl;
	for (Power j = 0; j <= _maxP; j++)
	  std::cout << this->operator()(i, j) << std::endl;
      }
    }
  protected:
    /**
     * Calculate index into internal storage for a given harmonic and
     * power.  @b NB No bounds check.
     *
     * @param n Harmonic
     * @param p Power
     *
     * @return Index into internal array
     */
    unsigned short index(Harmonic n, Power p) const
    {
#if 0
      if (p > _maxP) {
	std::cerr << "Index p=" << p << " > " << _maxP << std::endl;
	return 0xFFFF;
      }
      if (n > _maxN || -n < -_maxN) {
	std::cerr << "Index |n|=" << n << " > " << _maxN << std::endl;
	return 0xFFFF;
      }
#endif
      return ((n + _maxN) * (_maxP+1) + p);
    }
    Size           _maxN;         /**< Maximum harmonic order  */
    Size           _maxP;         /**< Maximum power of harmonic */
    bool           _useWeights;   /**< Wheter to use weights or not */
    ComplexVector  _q;            /**< Internal storage of Q vector */
  };
}
#endif
// Local Variables:
//  mode: C++
// End:
