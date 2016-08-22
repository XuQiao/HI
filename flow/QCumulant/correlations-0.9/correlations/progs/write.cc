/**
 * @file   correlations/progs/write.cc
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  An example of writing a data file.  
 * 
 * The program takes a number of options.  Do 
 * <pre class="shell">
 * ./write -h
 * </pre>
 * for information. 
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
#include <correlations/test/WriteData.hh>
#include <fstream>

/**
 * Show usage information
 *
 * @param prog Run name
 */
void
usage(const char* prog)
{
  using correlations::test::helpline;
  std::cout << "Usage: " << prog << " [OPTIONS]\n\n" << "Options:" << std::endl;

  helpline(std::cout, 'h', "", "This help", "");
  helpline(std::cout, 'o', "FILENAME", "Output file name", "data.dat");
  helpline(std::cout, 'e', "NEVENTS", "Number of events to write", "100");
  helpline(std::cout, 'm', "NPART", "Least number of particles/events", "800");
  helpline(std::cout, 'M', "NPART", "Largest number of particles/events",
           "1000");
}

/**
 * Entry point for program.
 *
 * Generate data file, and exit
 *
 * How many events to make
 * are set by the option @c -e, and the multiplicity range is
 * specified with options @c -m and @c -M.  The result is written to
 * the file given by the option @c -i, or @c data.dat if not specified
 *
 * @param argc Number of arguments
 * @param argv Vector of arguments
 *
 * @return 0 on success
 */
int
main(int argc,
     char** argv)
{
  unsigned short nEvents = 100;
  unsigned short minN = 800;
  unsigned short maxN = 1000;
  std::string output("data.dat");
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'h':
        usage(argv[0]);
        return 0;
      case 'e':
        nEvents = atoi(argv[++i]);
        break;
      case 'm':
        minN = atoi(argv[++i]);
        break;
      case 'M':
        maxN = atoi(argv[++i]);
        break;
      case 'o':
        output = argv[++i];
        break;
      default:
        std::cerr << argv[0] << ": Unknown option " << argv[i] << std::endl;
        return 1;
      }
    }
  }
  std::ofstream out(output.c_str());
  correlations::test::WriteData writer(minN, maxN);
  for (unsigned short ev = 0; ev < nEvents; ev++)
    writer.event(out, ev);
  out.close();
  return 0;
}
/*
 * EOF
 */
