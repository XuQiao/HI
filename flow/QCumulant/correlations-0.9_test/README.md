# Generic Multi-Particle Correlations


This code calculates the generic multi-particle correlator

	C{n} = <exp(i(sum_j^n h_j phi_j))>

where the vector `h=h_1,...,h_n` is the harmonics of each particle,
and `phi_j` is the observations of the azimuth angle.

See also 

    http://arxiv.org/abs/1312.3572

## Content:

* `cas` - Computer Algebra System files.
* `correlations` - C++ header library
* `data` - Some test data
* `doc` - Some code documentation files

## Building:

To check out, do (CERN users)

    svn co https://svn.cern.ch/reps/alicefmd/mcorrelations/trunk mcorrelations

For non-CERN users do

    svn co http://svnweb.cern.ch/guest/alicefmd/mcorrelations/trunk mcorrelations

Alternatively download from 

    http://www.nbi.dk/~cholm/mcorrelations/

and extract. 

To build, do

    cd mcorrelations
    make

To run tests, do

    make test
    make Test

(Test requires ROOT)

To make the documentation (requires Doxygen), do

	make doc

## Methods:

The code implements:

* [`correlations::NestedLoops`][nlb] - direct loop calculations.
* [`correlations::recursive::NestedLoops`][nlr] - loop calculations, using recursion.
* [`correlations::closed::FromQVector`][fqc] - Calculation from
  Q-vector using fixed expression (up to 8-particle correlator
  defined).
* [`correlations::recurrence::FromQVector`][fqr] - Calculation from
  Q-vector using combinatorics and recurrence
* [`correlations::recursive::FromQVector`][fqs] - Calculation from
  Q-vector using recursion only.

See also the example [`correlations/progs/analyze.cc`][ana] (and
[`correlations/progs/Analyze.C`][anr] - the ROOT version).

[nlb]: html/structcorrelations_1_1_nested_loops.html
[nlr]: html/structcorrelations_1_1recursive_1_1_nested_loops.html
[fqc]: html/structcorrelations_1_1closed_1_1_from_q_vector.html
[fqr]: html/structcorrelations_1_1recurrence_1_1_from_q_vector.html
[fqs]: html/structcorrelations_1_1recursive_1_1_from_q_vector.html
[ana]: html/analyze_8cc-example.html
[anr]: html/_analyze_8_c-example.html

## Copyright:


Copyright (c) 2013, Kristjan Gulbrandsen, Ante Bilandzic, Christian
Holm Christensen

## License:

GNU General Public License version 3 - see `COPYING`

