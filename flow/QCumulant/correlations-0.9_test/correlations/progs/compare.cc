/**
 * @file   correlations/progs/compare.cc
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 24 23:45:40 2013
 *
 * @brief  Compare outputs
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
#include <correlations/test/Comparer.hh>
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
  helpline(std::cout, 'a', "FILENAME", "First (A) file", "");
  helpline(std::cout, 'b', "FILENAME", "Second (B) file", "closed.dat");
}

/**
 * Entry point for program.
 *
 * Compare to output files of analyze
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
  std::string a("");
  std::string b("closed.dat");
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'h':
        usage(argv[0]);
        return 0;
      case 'a': a = argv[++i]; break;
      case 'b': b = argv[++i]; break;
      default:
        std::cerr << argv[0] << ": Unknown option " << argv[i] << std::endl;
        return 1;
      }
    }
  }
  std::ifstream af(a.c_str());
  std::ifstream bf(b.c_str());
  if (!af || !bf) {
      std::cerr << "Failed to open one or more of " << a << " and " << b << std::endl;
      return 1;
  }
  bool ret = correlations::test::Comparer::compare(std::cout, af, bf);
  return ret ? 0 : 1;
}
/*
 * EOF
 */
