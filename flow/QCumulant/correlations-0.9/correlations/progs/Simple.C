/*
 * This script can be AcLic compiled - i.e.,
 *
 *   Root> .L correlations/progs/Simple.C++
 *
 * If not, then the library headers will be.  
 */
/**
 * @file   Simple.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr  9 10:11:39 2014
 * 
 * @brief A simple ROOT example calculating n-particle correlations
 * for a single event with 11 particles.
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
#ifndef __CINT__
# include <correlations/Types.hh>
# include <correlations/Result.hh>
# include <correlations/QVector.hh>
# include <correlations/NestedLoops.hh>
# include <correlations/recursive/FromQVector.hh>
# include <correlations/recurrence/FromQVector.hh>
# include <correlations/closed/FromQVector.hh>
#else 
# include <vector>
#endif

/** 
 * A simple example.  This calculates up to the 8-particle correlator
 * for a single event with 11 particles using weights.  The
 * correlators are evaluated using 3 methods:
 *
 * - correlations::NestedLoops
 * - correlations::closed::FromQVector 
 * - correlations::recursive::FromQVector 
 *
 * The final results of each correlator, along with the time it took
 * to evaluate is printed at the end.
 * 
 */
void Simple()
{
  // --- Options -----------------------------------------------------
  // Whether to use weights 
  Bool_t useWeights = true;

  // Whether to calculte using nested loops 
  bool evaluateNestedLoops = true;

  // Maximim correlator order 
  size_t maxCorrelator = 6;

  // --- Data --------------------------------------------------------
  // Number of particles 
  const size_t nParticles      = 11;

  // Azimuthal angles:
  double angles[nParticles]  = { 6.28156, 1.02359, 1.77574, 
				 5.95144, 1.45554, 3.04718, 
				 6.01601, 4.67661, 3.39319, 
				 4.64926, 2.49226 };
  
  // Particle weights:
  double weights[nParticles] = { 1.14393, 0.90291, 0.93739, 
				 1.12879, 0.92277, 0.99567, 
				 1.13175, 1.07036, 1.01153, 
				 1.06911, 1.11211 };

  // --- Load code ---------------------------------------------------
#ifdef __CINT__
  if (maxCorrelator >= 7) 
    // Enable compilation of 7-particle correlator - slow
    gSystem->AddIncludePath("-DCORRELATIONS_CLOSED_ENABLE_U7");
  if (maxCorrelator >= 8) 
    // Enable compilation of 8-particle correlator  - slow
    gSystem->AddIncludePath("-DCORRELATIONS_CLOSED_ENABLE_U8");

  gROOT->LoadMacro("correlations/Result.hh+");
  gROOT->LoadMacro("correlations/QVector.hh+");
  gROOT->LoadMacro("correlations/Correlator.hh+");
  gROOT->LoadMacro("correlations/NestedLoops.hh+");
  // Work around bug in ROOT 
  gSystem->CopyFile("correlations/recursive/FromQVector.hh",
		    "correlations/recursive/_FromQVector.hh");
  gROOT->LoadMacro("correlations/recursive/_FromQVector.hh+");
  gROOT->LoadMacro("correlations/closed/FromQVector.hh+");
  gROOT->LoadMacro("correlations/test/Stopwatch.hh+");
  // Remove work-around temporary 
  gSystem->Unlink("correlations/recursive/_FromQVector.hh");
#endif


  // --- Derived data ------------------------------------------------
  // Our harmonics 
  correlations::Size maxH = 8;
  correlations::HarmonicVector h(maxH);
  h[0] = 1; 
  h[1] = 3; 
  h[2] = 5; 
  h[3] = 0; 
  h[4] = -2; 
  h[4] = -4; 
  h[6] = -1; 
  h[7] = -6;

  // Our Q-vector 
  correlations::QVector q(h,useWeights);
  
  // Our correlators 
  correlations::closed::FromQVector closed(q);
  correlations::recursive::FromQVector recursive(q);
  correlations::RealVector ap(nParticles);
  correlations::RealVector aw(nParticles);
  correlations::NestedLoops nested(ap,aw);
  correlations::Correlator* correlators[] = { 
    &closed, &recursive, (evaluateNestedLoops ? &nested : 0), 0 };

  // --- Particle loop -----------------------------------------------
  // Calculate Q-vector 
  for (size_t iParticle = 0; iParticle < nParticles; iParticle++) { 
    ap[iParticle] = angles[iParticle];
    aw[iParticle] = weights[iParticle];
    q.fill(angles[iParticle], weights[iParticle]);
  }

  // --- Calculation loop --------------------------------------------
  correlations::test::Stopwatch* s = correlations::test::Stopwatch::create();
  for (correlations::Size i = 2; i <= maxCorrelator; i++) {
    printf("%d-particle correlator:\n", i);
    
    correlations::Correlator** p = correlators;
    while (*p) { 
      s->start(true);
      correlations::Correlator* c  = (*p);
      correlations::Result      r  = c->calculate(i, h);
      correlations::Complex     rc = r.eval();
      s->stop();
      printf(" %-22s: %12g + %12gi (%12gs)\n", 
	     c->name(), rc.real(), rc.imag(), s->seconds());
      p++;
    }
  }
}
// EOF

    
    
