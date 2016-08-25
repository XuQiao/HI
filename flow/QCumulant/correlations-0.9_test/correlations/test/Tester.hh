#ifndef CORRELATIONS_TEST_TESTER_HH
#define CORRELATIONS_TEST_TESTER_HH
/**
 * @file   correlations/test/Tester.hh
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
#include <correlations/recursive/FromQVector.hh>
#include <correlations/recursive/NestedLoops.hh>
#include <correlations/recurrence/FromQVector.hh>
#include <correlations/closed/FromQVector.hh>
#include <correlations/test/Random.hh>
#include <correlations/test/ReadData.hh>
#include <correlations/test/Printer.hh>
#include <correlations/test/Stopwatch.hh>

namespace correlations
{
  namespace test
  {
    /**
     * A class to test the code.
     *
     * The input is assumed to have been generated by
     * correlations::test::WriteData.  The result is written to a
     * stream - typically the screen.  The results can also be saved
     * into a stream for later comparison.
     *
     */
    struct Tester
    {
      /**
       * Modes
       */
      enum EMode
      {
        CLOSED, RECURRENCE, RECURSIVE
      };
      /**
       * @param s Input string
       *
       * @return Mode id
       */
      static EMode
      str2mode(const std::string& s)
      {
        if (s == "CLOSED")
          return CLOSED;
        else if (s == "RECURRENCE")
          return RECURRENCE;
        else if (s == "RECURSIVE")
          return RECURSIVE;
        std::cerr << "Unknown mode: " << s << " assuming CLOSED" << std::endl;
        return CLOSED;
      }
      /**
       * Constructor
       *
       * @param input    Input stream
       * @param mode     What algorithms to use
       * @param maxN     Max # of particles to correlate
       * @param doNested Whether to run nested loop code
       * @param verbose  Whether to be verbose
       */
      Tester(std::istream& input, EMode mode = CLOSED, Size maxN = 8,
          bool doNested = false, bool verbose = false) :
          _h(maxN), _phis(), _weights(), _r(input), _q(0, 0, true), _c(0), _n(
              0), _rC(0), _rN(0), _s(0), _tC(0), _tN(0), _e(0), _v(verbose)
      {
        if (_v)
          std::cout << "Harmonics:" << std::flush;
        Size sum = 0;
        Random::seed(54321);
        for (Size i = 0; i < maxN; i++) {
          _h[i] = Random::asHarmonic(-6, 6);
          if (_v)
            std::cout << " " << _h[i];
          sum += std::abs(_h[i]);
        }
        if (_v)
          std::cout << " -> " << sum << std::endl;
        _q.resize(_h);
        _rC.resize(maxN - 1);
        _rN.resize(doNested ? maxN - 1 : 0);
        _tC.resize(maxN - 1);
        _tN.resize(doNested ? maxN - 1 : 0);
        _s = Stopwatch::create();
        switch (mode)
          {
        case RECURSIVE:
          _c = new correlations::recursive::FromQVector(_q);
          break;
        case RECURRENCE:
          _c = new correlations::recurrence::FromQVector(_q);
          break;
        case CLOSED:
          _c = new correlations::closed::FromQVector(_q);
          break;
          }
        if (doNested)
          switch (mode)
            {
          case RECURSIVE:
          case RECURRENCE:
            _n = new correlations::recursive::NestedLoops(_phis, _weights,
                true);
            break;
          case CLOSED:
            _n = new correlations::NestedLoops(_phis, _weights, true);
            break;
            }
        if (_v)
          {
            std::cout << "Cumulant correlator: " << _c->name() << std::endl;
            if (_n)
              std::cout << "Nested loop correlator: " << _n->name()
                  << std::endl;
          }
      }
      /**
       * Destructor
       */
      virtual
      ~Tester()
      {
      }
      /**
       * Make a single event
       *
       */
      bool
      event()
      {
        _q.reset();
        if (!_r.event(_q, _phis, _weights))
          return false;

        _e++;
        if (_v)
          std::cout << "Event # " << std::setw(4) << _e << ": " << std::setw(4)
              << _phis.size() << " particles " << std::flush;
        for (Size i = 0; i < _rC.size(); i++)
          {
            Size n = i + 2;
            if (_v)
              std::cout << (i == 0 ? "" : "..") << n << std::flush;
            _s->start(true);
            _rC[i] += _c->calculate(n, _h);
            _tC[i] += _s->stop();
            // _s->Print();
            // std::cout << "Calculated QC" << std::endl;

            if (_n)
              {
                if (_v)
                  std::cout << '+' << std::flush;
                _s->start(true);
                _rN[i] += _n->calculate(n, _h);
                _tN[i] += _s->stop();
              }
            // std::cout << "Calculated nested loops" << std::endl;
          }
        if (_v)
          std::cout << " done" << std::endl;

        return true;
      }
      /**
       * Do calculations at the end
       *
       * @param out  Where to print the results
       */
      void
      end(std::ostream& out)
      {
        Printer::title(out, "From Q-vector", "From loops", "T_Q", "T_loops");
        for (Size i = 0; i < _rC.size(); i++)
          {
            Complex rc = _rC[i].eval();
            Complex rn = _n ? _rN[i].eval() : Complex(0, 0);
            Real tc = _tC[i] / _e;
            Real tn = _n ? _tN[i] / _e : -1;
            Printer::result(out, 2 + i, rc, rn, tc, tn);
          }
      }
      /**
       * Save results to file
       *
       * @param out Output file
       */
      void
      save(std::ostream& out)
      {
        size_t savePrec = out.precision();
        out.precision(16);
        out << "# A total of " << _rC.size() << " cumulants\n"
            << "# Harmonics:\n"
            << _h.size();
        for (Size i = 0; i < _h.size(); i++) out << std::setw(5) << _h[i];
        out << "\n# Order   QC    NL     t_QC     t_NL" << std::endl;
        for (Size i = 0; i < _rC.size(); i++)
          {
            Complex rc = _rC[i].eval();
            Complex rn = _n ? _rN[i].eval() : Complex(0, 0);
            Real tc = _tC[i] / _e;
            Real tn = _n ? _tN[i] / _e : -1;
            out << i + 2 << "\t" << rc << "\t" << rn << "\t" << tc << "\t" << tn
                << std::endl;
          }
        out << "# EOF" << std::endl;
        out.precision(savePrec);
      }
    protected:
      Tester(const Tester&);
      Tester& operator=(const Tester&);
      /** Harmonics to use */
      HarmonicVector _h;
      /** Phi cache */
      RealVector _phis;
      /** Weight cache */
      RealVector _weights;
      /** My generator */
      ReadData _r;
      /** Our Q vector */
      QVector _q;
      /** Correlator that uses cumulants */
      FromQVector* _c;
      /** Correlator that uses nested loops */
      NestedLoops* _n;
      /** Cumulant results */
      ResultVector _rC;
      /** Nested loop results */
      ResultVector _rN;
      /** Stop watch */
      Stopwatch* _s;
      /** Cumulant timing */
      RealVector _tC;
      /** Nested loop timing */
      RealVector _tN;
      /** Counter of events */
      Size _e;
      /** Be verbose */
      bool _v;
    };
  }
}

#endif /* TESTER_HH_ */
