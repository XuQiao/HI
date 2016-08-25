#ifndef CORRELATIONS_TEST_STOPWATCH_H
#define CORRELATIONS_TEST_STOPWATCH_H
/**
 * @file   correlations/test/Stopwatch.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Sun Oct 27 21:58:37 2013
 *
 * @brief  Stop-watches
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
#include <ctime>
#if defined(__unix__) || defined(__unix) || defined(unix)
# include <unistd.h>
# if _POSIX_C_SOURCE >= 199309L
#  ifdef _POSIX_TIMERS
#   define HAVE_POSIX_TIMERS
#  endif 
# endif
#endif

namespace correlations {
  namespace test {
    /**
     * Base class for stop-watches
     *
     * @code
     * correlations::test::Stopwatch* s = correlations::test::Stopwatch::Create();
     * @endcode
     * returns the best possible stop-watch
     *
     * @headerfile "" <correlations/test/Stopwatch.hh>
     */
    struct Stopwatch
    {
      /**
       * Destructor
       */
      virtual ~Stopwatch() {}
      /**
       * Create a stop watch
       */
      static Stopwatch* create();
      /**
       * Test if the stop watch is running
       *
       * @return true if running
       */
      bool isRunning() const { return _state == kRunning; }
      /**
       * Test if the stop watch is stopped
       *
       * @return true if non-running
       */
      bool isStopped() const { return _state == kStopped; }
      /**
       * Reset the stop-watch
       */
      void reset()
      {
	_state = kStopped;
	doReset();
      }
      /**
       * Start the stop-watch
       *
       * @param rst If true, reset stop-watch before starting
       */
      void start(bool rst)
      {
	if (rst) reset();
	if (isRunning()) return;
	_state = kRunning;
	doStart();
      }
      /**
       * Stop the stop-watch.
       *
       * @return accumulated time in seconds
       */
      double stop()
      {
	if (isStopped()) return seconds();
	doStop();
	_state      =  kStopped;
	return seconds();
      }
      /**
       * Get the accumulated elapsed time in seconds
       *
       * @return Accumulated elapsed time in seconds
       */
      double seconds() const { return doSeconds(); }
    protected:
      /** The states */
      enum State_t {
	kRunning,
	kStopped
      } _state;
      /**
       * Constructor
       */
      Stopwatch() : _state(kStopped) {}
      /**
       * @{
       * @name Implementation interface
       */
      virtual void doStart() = 0;
      virtual void doStop() = 0;
      virtual void doReset() = 0;
      virtual double doSeconds() const = 0;
      /* @} */
    };

    /**
     * Name space for ISO code
     *
     */
    namespace ISO {
      /**
       * A stop watch using the standard clock() function.
       *
       * @note @c clock() is part of the ISO C standard, while the @c
       * clock_gettime() function is POSIX, so we opt for @c clock()
       * for portability reasons.
       *
       * This has a precision of 10's of miliseconds
       *
       * @headerfile "" <correlations/test/Stopwatch.hh>
       */
      struct Stopwatch : public correlations::test::Stopwatch
      {
	/**
	 * Constructor
	 */
	Stopwatch() : correlations::test::Stopwatch(), _start(0), _sum(0) {}
      protected:
	/**
	 * Reset the stop-watch
	 */
	void doReset() { _start = 0; _sum = 0;}
	/**
	 * Start the stop-watch
	 */
	void doStart() { _start = clock(); }
	/**
	 * Stop the stop-watch.
	 */
	void doStop()
	{
	  _sum += (clock() - _start);
	}
	double doSeconds() const
	{
	  return double(_sum) / (CLOCKS_PER_SEC);
	}
	/** When the timer was last started */
	clock_t _start;
	clock_t _sum;
      };
    }

#ifdef HAVE_POSIX_TIMERS
    /**
     * Name space for POSIX code
     *
     */
    namespace POSIX
    {
      /**
       * A stop watch using the POSIX clock_gettime() function.  This
       * has a precision of nano-seconds
       *
       * @headerfile "" <correlations/test/Stopwatch.hh>
       */
      struct Stopwatch : public correlations::test::Stopwatch
      {
	/**
	 * Constructor
	 */
	Stopwatch() : correlations::test::Stopwatch(), _start(), _sum() {}
      protected:
	/**
	 * Reset the stop-watch
	 */
	void doReset()
	{
	  _sum.tv_sec  = _start.tv_sec  = 0;
	  _sum.tv_nsec = _start.tv_nsec = 0;
	}
	/**
	 * Start the stop-watch
	 */
	void doStart()
	{
	  clock_gettime(_POSIX_CPUTIME, &_start);
	}
	/**
	 * Stop the stop-watch.
	 */
	void doStop()
	{
	  struct timespec ts;
	  clock_gettime(_POSIX_CPUTIME, &ts);
	  _sum.tv_sec   += (ts.tv_sec  - _start.tv_sec);
	  _sum.tv_nsec  += (ts.tv_nsec - _start.tv_nsec);
	}
	/**
	 * Get the accumulated number of seconds
	 *
	 * @return accumulated number of seconds
	 */
	double doSeconds() const
	{
	  return _sum.tv_sec + _sum.tv_nsec / 1e9;
	}
	/** When the timer was last started */
	struct timespec _start;
	struct timespec _sum;
      };
    }
#endif
    Stopwatch* Stopwatch::create()
    {
#ifdef HAVE_POSIX_TIMERS
      return new POSIX::Stopwatch();
#else
      std::cerr << "Warning, using poor ISO stopwatch" << std::endl;
      return new ISO::Stopwatch();
#endif
    }
  }
}
#endif
// Local Variables:
//  mode: C++
// End:
