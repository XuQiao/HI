#ifndef CORRELATIONS_TEST_READDATA_H
#define CORRELATIONS_TEST_READDATA_H
/**
 * @file   correlations/test/ReadData.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  Code to read in data
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
#include <correlations/QVector.hh>
#include <fstream>
#include <sstream>
#include <string>

namespace correlations {
  namespace test {
    /**
     * A Read data from file.
     *
     * @code
     * std::ifstream in("data.dat");
     * ReadData   reader(in);
     * QVector    q;
     * RealVector phi;
     * RealVector weight;
     * bool       more = true;
     *
     * do {
     *   more = reader.Event(q, phi, weight);
     * } while(more);
     * @endcode
     *
     * @headerfile "" <correlations/test/ReadData.hh>
     */
    struct ReadData
    {
      /**
       * Constructor
       *
       * @param input Stream to read from
       */
      ReadData(std::istream& input)
	: _input(input), _phiR(0)
      {
      }
      virtual ~ReadData() {}
      /**
       * Create one event
       *
       * @param q        Q-Vector to fill
       * @param phis     @f$ \phi@f$ vector to fill
       * @param weights  Weight vector to fill
       *
       * @return true if more events to process
       */
      virtual bool event(QVector& q,
			 RealVector& phis,
			 RealVector& weights)
      {
	if (_input.eof()) return false;

	std::string l;
	do {
	  std::getline(_input, l);
	  if (l.empty()) return false;
	  if (l[0] == '#') continue;
	  break;
	} while (true);

	// std::cout << "Event header " << l << std::endl;
	std::stringstream s1(l);
	int id = 0;
	s1 >> id;
	if (id != -1)
	  // Not a header
	  return false;

	Size mult = 0;
	Size ev   = 0;
	s1 >> mult >> ev >> _phiR;

	if (mult <= 0)
	  // Bad multiplicity
	  return false;

	phis.resize(mult);
	weights.resize(mult);
	do {
	  std::getline(_input, l);
	  if (l.empty()) return false;
	  if (l[0] == '#') continue;

	  std::stringstream s2(l);
	  s2 >> id;
	  if (id < 0) {
	    // std::cout << "Event trailer: " << l << std::endl;
	    break;
	  }

	  Real phi=0, weight=0;
	  s2 >> phi >> weight;
	  phis[id]    = phi;
	  weights[id] = weight;

	  q.fill(phi, weight);
	} while (true);
	if (id != -2)
	  // Not a trailer
	  return false;

	return true;
      }
      Real phiR() const {return _phiR; }
    protected:
      std::istream& _input;
      Real _phiR;
    };
  }
}
#endif
// Local Variables:
//  mode: C++
// End:
