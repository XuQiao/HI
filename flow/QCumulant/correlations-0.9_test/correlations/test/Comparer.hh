#ifndef CORRELATIONS_TEST_Comparer_HH
#define CORRELATIONS_TEST_Comparer_HH
/**
 * @file   correlations/test/Comparer.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  An example
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
#include <correlations/test/Printer.hh>

namespace correlations {
  namespace test {
    /** 
     * Compare results from two streams. 
     *
     * The results are assumed to come from
     * correlations::test::Tester.  This class will then show the
     * difference between the results as well as the difference in
     * timing.
     */
    struct Comparer
    {
      /**
       * Compare two results files
       *
       * @param out Output stream
       * @param in1 First input file
       * @param in2 Second input file
       *
       * @return true if results match
       */
      static bool
      compare(std::ostream& out,
              std::istream& in1,
              std::istream& in2)
      {
        HarmonicVector h1;
        HarmonicVector h2;
        ComplexVector r1;
        ComplexVector r2;
        RealVector t1;
        RealVector t2;

        if (!readFile(in1, h1, r1, t1)) {
          std::cerr << "Failed to read file 1" << std::endl;
          return false;
        }
        if (!readFile(in2, h2, r2, t2)) {
          std::cerr << "Failed to read file 2" << std::endl;
          return false;
        }
        Size maxH = std::max(h1.size(), h2.size());
        if (maxH <= 1) {
          std::cerr << "Nothing to compare" << std::endl;
          return false;
        }
        maxH--;
        if (h1.size() != h2.size()) std::cerr
            << "Inconsistent number of harmonics, only using " << maxH
            << std::endl;

        bool ret = true;
        for (Size i = 0; i <= maxH; i++) {
          if (h1[i] != h2[i]) {
            std::cerr << "Harmonic # " << i + 2 << " inconsistent: " << h1[i]
                      << " != " << h2[i] << std::endl;
            ret = false;
          }
        }

        Printer::title(out, "A", "B", "(T_B-T_A)/T_B");
        for (Size i = 0; i < maxH; i++) {
          Real dt = (t2[i] - t1[i]) / t1[i];
          if (Printer::result(out, 2 + i, r1[i], r2[i], dt)) ret = false;
        }
        return ret;
      }
      /**
       * Read in a file
       *
       * @param in File to read
       * @param h  On return, the harmonics
       * @param r  On return, the results read
       * @param t  On return. the timing read
       *
       * @return true if the file was read
       */
      static bool
      readFile(std::istream& in,
               HarmonicVector& h,
               ComplexVector& r,
               RealVector& t)
      {

        while (!in.eof()) {
          std::string l;
          std::getline(in, l, '\n');
          if (in.bad()) {
            std::cerr << "Error while reading line from comparison file"
                      << std::endl;
            return false;
          }
          if (l.empty()) continue;
          if (l[0] == '#') continue;

          std::stringstream sf(l);
          if (h.size() <= 0) {
            Size nh;
            sf >> nh;
            if (nh <= 0) {
              std::cerr << "Invalid number of harmonics: " << nh << std::endl;
              return false;
            }
            h.resize(nh);
            for (Size j = 0; j < nh; j++)
              sf >> h[j];
            r.resize(h.size());
            t.resize(h.size());
            continue;
          }
          Size i;
          Complex rc;
          Complex rn;
          Real tc;
          Real tn;
          sf >> i >> rc >> rn >> tc >> tn;
          if (sf.bad()) {
            std::cerr << "Error while extracting from line " << sf.str()
                      << std::endl;
            return false;
          }

          if (i < 2) continue;
          Size n = i - 2;
          if (n >= h.size()) continue;

          r[n] = rc;
          t[n] = tc;
        }
        return true;
      }
    };
  }
}

#endif /* COMPARE_HH_ */
//
// EOF
//
