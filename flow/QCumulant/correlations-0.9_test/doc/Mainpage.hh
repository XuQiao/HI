/**
 * @file   Mainpage.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Nov 13 11:03:21 2013
 *
 * @brief  Documentation file
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
#ifndef DOXYGEN_INPUT
#error Not for compilation
#endif
/**
 * @mainpage Calculate n-particle Correlations
 *
 * @tableofcontents
 *
 * @section intro Introduction
 *
 * This code can calculate the _n_-particle correlator
 *
 * @f[
 *  C\{N\}=
 *  \left\langle\exp\left[i\left(\sum_j^N h_j\phi_j\right)\right]\right\rangle
 * @f]
 *
 * for any @f$ N@f$ and harmonic orders @f$ h_i@f$ - limited by
 * computing time and available memory.
 *
 * See also 
 * <a href="http://arxiv.org/abs/1312.3572">arXiv:1312.3572 [nucl-ex]</a>.
 *
 * @section downinst Download and Installation
 *
 * @subsection checkout Checkout from Subversion
 *
 * The code is managed by Subversion.  The repository can be found at
 * http://svnweb.cern.ch/world/wsvn/alicefmd/mcorrelations/
 *
 * The code can be checked out by (may require CERN credentials)
 *
 * <pre class="shell">
 * svn co https://svn.cern.ch/reps/alicefmd/mcorrelations/trunk mcorrelations</pre>
 *
 * If the above does not work, try 
 *
 * <pre class="shell">
 * svn co http://svnweb.cern.ch/guest/alicefmd/mcorrelations/trunk mcorrelations</pre>
 *
 *
 * @subsection download Download
 *
 * Alternatively, the code can be downloaded from
 * http://www.nbi.dk/~cholm/mcorrelations/correlations-0.9.tar.gz
 *
 * @subsection building Building the package
 *
 * All code of the project is contained in declaration (header) files.
 * That is, all code is 'inline'.  That means that there is no need to
 * compile a library or the like first.
 *
 * @subsection usage_other Using in Another Project 
 * 
 * To use the code in some other project located at say
 * @c /some/directory do
 *
 * <pre class="shell">make install PREFIX=/some/directory</pre>
 *
 * and then make sure that @c /some/directory is in the include path
 * of the compiler (e.g., <tt>-I/some/directory</tt>).  One can also
 * copy the sub-directory <tt>correlations</tt> to the project in
 * question and use it directly from there (please observe the license
 * conditions).
 *
 * @note Concerning compatibility.  The code it self is portable as it
 * is, and only uses ISO C++98/C99.  However, the code and in
 * particular the build system captured in the top-level Makefile have
 * only been tested on GNU/Linux with relatively new versions of GCC.
 * To build the tests on other platforms it may require a bit of
 * tweaking on the users behalf.
 *
 * - 
 *
 * @note When using the code with <a href="root.cern.ch">ROOT</a> it
 * is important that the code is (at least) <a
 * href="http://root.cern.ch/drupal/content/compiling-your-code-also-known-aclic">ACLiC</a>
 * compiled.  Suppose we have the ROOT script <tt>Test.C</tt> @include
 * Test.C and the code of this package is in <tt>/some/directory</tt>,
 * then we need to do
 * <pre class="shell">
 * Root> gSystem->AddIncludePath("-I/some/directory");
 * Root> .L Test.C+
 * Root> Test() </pre>
 *
 * @subsection tests Examples and Tests
 *
 * The package comes with a number of examples and tests
 *
 * The test code lives in the namespace correlations::test:
 *
 * - correlations::test::Random generates random numbers
 * - correlations::test::Distribution allows one to sample a
 *   @f$\phi@f$ distribution.
 * - correlations::test::Weights generates weights
 * - correlations::test::WriteData writes a data file
 * - correlations::test::ReadData reads in a data file
 * - correlations::test::Stopwatch times execution of code blocks
 * - correlations::test::Tester tests the correlation code using an
 *   existing data file
 * - correlations::test::Printer prints results to a stream (a file or
 *   the screen).
 * - correlations::test::Comparer compares the results of two runs of
 *   correlations::test::Tester
 *
 * There are some predefined examples that uses these classes.  The
 * code is in the sub-directory correlations/prog:
 *
 * - <a href="write_8cc-example.html">write.cc</a> generates a data file
 * - <a href="analyze_8cc-example.html">analyze.cc</a> analyses the
 *   data file generated for correlations
 * - <a href="compare_8cc-example.html">compare.cc</a> compares the
 *   results of two different runs of analyze.cc
 *
 * To build and run the tests, do
 *
 * <pre class="shell">
 * make test
 * </pre>
 *
 * @note By default, only up to 6-particle correlators are defined for
 * the so-called 'fixed-form' (see @ref algo), but 7 and 8-particle
 * correlators can be enabled by doing
 * <pre class="shell">
 * make test USE7=1 USE=8
 * </pre>
 * Do note, that compilation will take a very long time in this case.
 *
 * All of these programs allow for a set of options.  To see which
 * options are accepted, do
 * <pre class="shell">
 * ./prog -h
 * </pre>
 * where <tt>prog</tt> is one of <tt>write</tt>, <tt>analyze</tt>, or
 * <tt>compare</tt>.
 *
 * The ROOT script <a href="_simple_8_c-example.html">Simple.C</a> is
 * a simple example that calculates correlations for a single event.
 *
 * In addition, there are 3 ROOT-based scripts/programs to illustrate
 * the use of the code in ROOT.  These are parallel to the 3 above,
 * and are:
 *
 * - <a href="_write_8_c-example.html">Write.C</a> generates a data file
 * - <a href="_analyze_8_c-example.html">Analyze.C</a> analyses the
 *   data file generated for correlations
 * - <a href="_compare_8_c-example.html">Compare.C</a> compares the
 *   results of two different runs of analyze.cc
 *
 * To build and run these examples do
 * <pre class="shell">
 * make Test
 * </pre>
 * (notice the change of case)
 *
 * @note By default, only up to 6-particle correlators are defined for
 * the so-called 'fixed-form' (see @ref algo), but 7 and 8-particle
 * correlators can be enabled by doing
 * <pre class="shell">
 * make Test USE7=1 USE=8
 * </pre>
 * Do note, that compilation will take a very long time in this case.
 *
 * The scripts can also be executed in the interpretive environment -
 * but should always be <a
 * href="http://root.cern.ch/drupal/content/compiling-your-code-also-known-aclic">ACLiC</a>:
 *
 * <pre class="shell">
 * Root> gSystem->AddIncludePath("-I.");
 * Root> .x correlations/test/Write.C+(nEvents,minN,maxN)
 * Root> .x correlations/test/Analyze.C+("CLOSED",maxH,doNestedLoops)
 * Root> .x correlations/test/Analyze.C+("RECURRENCE",maxH,doNestedLoops)
 * Root> .x correlations/test/Analyze.C+("RECURSIVE",maxH,doNestedLoops)
 * Root> .x correlations/test/Compare.C+("CLOSED","RECURRENCE")
 * Root> .x correlations/test/Compare.C+("CLOSED","RECURSIVE")
 * </pre>
 *
 * @section algo Algorithms
 *
 * The correlation is calculated using the concrete sub-classes of
 * correlations::Correlator:
 *
 * - correlations::NestedLoops uses nested loops (or rather, recursive
 *   function calls.  This relies on storing all the observations of
 *   @f$ \phi@f$ and (possibly) the corresponding weights @f$ w@f$.
 *   This way of calculating the correlation scales like @f$ O(M^N)@f$
 *   where @f$ M@f$ is the event multiplicity and @f$ N@f$ is the
 *   correlation order to calculate.  It is indeed rather slow.
 *
 * - correlations::FromQVector used a @f$ Q@f$-vector to calculate the
 *   correlations.  Calculating the @f$ Q-vector@f$ is an order @f$
 *   O(N M)@f$ operation, while calculating the correlation is only
 *   @f$ O(N)@f$.  Still - it is much faster generally, then the
 *   nested loop method.
 *
 * Both algorithms are interfaced in the same way through the member
 * function correlations::Correlator::calculate.
 *
 * The member function correlations::Correlator::calculate returns a
 * correlations::Result object, which contains the sum and the sum of
 * the weights.  The sum and the sum of weights can be accumulated in
 * an object of this time, and the final event average result can be
 * extracted using correlations::Result::eval.
 *
 * @subsection recursion Recursive vs. closed-form
 *
 * The class correlations::FromQVector is again an abstract base class
 * (cannot make objects of that type), and two different
 * implementations are provided.
 *
 * - correlations::recurrence::FromQVector which used recursion to
 *   simplify the expression, and has the added feature that it can
 *   calculate any order correlator - provided enough computing time
 *   and memory.
 * - correlations::closed::FromQVector which uses closed-form
 *   expression for the correlators.  This means that we are limited
 *   to the order for which we have generated these closed forms.
 *   Note, one can use the program print.cc to generate this, should
 *   it be needed.
 *
 * Similarly to correlations::recursive::FromQVector, there is a
 * sub-class of correlations::NestedLoops -
 * correlations::recursive::NestedLoops which does the nested loops
 * using recursion.
 *
 * @note The recursive algorithms could probably be implemented using
 * template recursion to let the compiler (and not run-time) deal with
 * branching, etc., which could speed up this tremendously.
 *
 * @section qvec Q-Vector
 *
 * The class correlations::QVector is the container of the per-event
 * @f$Q@f$-vector
 *
 * @f[
 *  Q(h,p) = \sum_j w_j^p e^{ih\phi_j}
 * @f]

 * where @f$ \phi_j@f$ is the angle of the observation, and @f$ w_j@f$
 * it's weight. Observations are added to the @f$ Q@f$-vector using
 * the member function correlations::QVector::fill.
 *
 * @section impl Implementation notes.
 *
 * - All code lives in the namespace correlations. This to minimize
 *   the pollution of users code space.
 * - the file correlations/Types.hh defines a number of basic types,
 *   like correlations::Real and correlations::Complex and the rest of
 *   the code uses these type definitions rigeroursly.  This allows
 *   one to easily change the underlying types should the need be -
 *   e.g., changing correlations::Real to be an alias for long double
 *   rather than double would allow for greater precision.
 * - The methods are implemented in plain C++.
 * - Testing code is available in correlations::test
 *
 * @section liclink Authors, Copyright, and License
 *
 * Copyright &copy; 2013 Kristjan Gulbrandsen, Ante Bilandzic,
 * Christian Holm Christensen.
 *
 * The code was developed by
 *
 * - Kristjan Gulbrandsen <a href="mailto:gulbrand@nbi.dk">gulbrand AT
 *   nbi DOT dk</a>
 * - Ante Bilandzic <a href="mailto:bilandz@nbi.dk">bilandz AT nbi DOT
 *   dk</a>
 * - Christian Holm Christensen <a href="mailto:cholm@nbi.dk">cholm AT
 *   nbi DOT dk</a>
 *
 * The code is distributed under the GNU General Public License
 * version 3.  See the page @ref license
 *
 * @example Simple.C A simple example 
 *
 * @example write.cc A simple program that creates a data file using
 * correlations::test::Distribution correlations::test::Weights and
 * correlations::test::WriteData
 *
 * @example analyze.cc A simple program that runs the algorithms over
 * a data file.  The test is steered by correlations::test::Tester and
 * the data read by correlations::test::ReadData.
 *
 * @example compare.cc A simple program to compare the results of
 * running analyze.cc twice. Uses correlations::test::Comparer
 *
 * @example Write.C A ROOT-based script to generate a data file
 *
 * @example Analyze.C A ROOT-based script to illustrate the use of
 * the code in ROOT
 *
 * @example Compare.C A ROOT-based script to compare the results of
 * running Analyze.C twice.
 *
 * @example print.cc A simple program that dumps the expressions for
 * the correlations using @f$ Q@f$-vector input and recursion.
 *
 */

/*
 * EOF
 */
