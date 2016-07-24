#ifndef CORRELATIONS_RESULT_HH
#define CORRELATIONS_RESULT_HH
/**
 * @file   correlations/Result.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  Results of calculations
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
#include <iostream>

namespace correlations {
  //____________________________________________________________________
  /**
   * Structure to hold and accummulate the results in
   *
   @code
   correlations::QVector q;
   correlations::Result  r;
   correlations::Correlator c = ...;
   correlations::HarmonicVector h(n);
   
   while (moreEvents) {
   
     while (moreObservations) {
       correlations::Real phi    = NextObservation();
       correlations::Real Weight = GetWeight(phi);
   
       q.fill(phi, weight);
     }
   
     r += c.calculate(h);
   }
   std::cout << r.eval() << std::endl;
   @endcode
   * @headerfile ""  <correlations/Result.hh>
   */
  struct Result
  {
    Result() : _sum(0,0), _weights(0) {}
    /**
     * Constructor
     *
     * @param c Complex number
     * @param w Weigth
     */
    Result(const Complex& c, Real w)
      : _sum(), _weights(w)
    {
    	_sum = c;
    }
    /**
     * Assignment plus operator
     *
     * @param r Result to add to this
     *
     * @return Reference to this
     */
    Result& operator+=(const Result& r)
    {
      _sum  += r._sum;
      _weights += r._weights;
      return *this;
    }
    Result& operator=(const Complex& c)
    {
      _sum = c;
      return *this;
    }
    /**
     * Evaluate the result
     *
     * @return The scaled result
     */
    Complex eval() const
    {
      return (_weights != 0 ? _sum / _weights : 0);
    }
    void print() const
    {
      std::cout << _sum << "\t" << _weights << std::endl;
    }
 // protected:
    /** Sum of arguments */
    Complex _sum;
    /** Sum of Weights */
    Real  _weights;
  };
  /** Vector of results */
  typedef std::vector<Result> ResultVector;
}
#endif
// Local Variables:
//  mode: C++
// End:
