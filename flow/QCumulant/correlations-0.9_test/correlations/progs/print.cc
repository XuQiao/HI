/**
 * @file   correlations/progs/print.cc
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Mon Oct 28 10:39:32 2013
 *
 * @brief  Code to generate recursive correlator functions
 *
 * Compile like
 *
 * @verbatim
 *   g++ print.cc -o print
 * @verbatim
 *
 * Run like
 *
 * @verbatim
 *   ./print -n 8 -t latex
 *   ./print -n 5 -t C++
 *   ./print -n 12 -t maxima
 * @endverbatim
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
#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <cmath>
/** Type of real variables */
typedef double Real;
/** Type of powers */
typedef unsigned short int Power;
/** type of sizes */
typedef unsigned short int Size;
/** Vector of sizes */
typedef std::vector<Size> SizeVector;

//____________________________________________________________________
/**
 * Return the next combination of @a k element from the range
 * [@a first, @a last).
 *
 * @param first  Start of the range
 * @param k      Number of elements to get
 * @param last   End of range
 *
 * @return true if there are more (unique) combinations
 */
bool next_combination(const SizeVector::iterator first,
		      SizeVector::iterator       k,
		      const SizeVector::iterator last)
{
  /* Credits: Mark Nelson http://marknelson.us */
  if ((first == last) || (first == k) || (last == k)) return false;
  SizeVector::iterator i1 = first;
  SizeVector::iterator i2 = last;
  ++i1;
  if (last == i1) return false;

  i1 = last;
  --i1;
  i1 = k;
  --i2;

  while (first != i1) {
    if (*--i1 < *i2) {
      SizeVector::iterator j = k;
      while (!(*i1 < *j)) ++j;
      std::iter_swap(i1,j);
      ++i1;
      ++j;
      i2 = k;
      std::rotate(i1,j,last);
      while (last != j) {
	++j;
	++i2;
      }
      std::rotate(k,i2,last);
      return true;
    }
  }
  std::rotate(first,k,last);
  return false;
}

//____________________________________________________________________
/**
 * Printer base class
 */
struct Printer
{
  Printer(std::ostream& out) : mOut(out) {}
  virtual ~Printer() {}
  /**
   * Print out stuff before doing any expansions
   *
   * @return Reference to this
   */
  virtual Printer& start() { return *this; }
  /**
   * Print out stuff after expansion
   *
   * @param Size Number of declarations written
   *
   * @return Reference to this
   */
  virtual Printer& stop(Size) { return *this; }
  /**
   * Start of a function
   *
   * @param begin Iterator pointing to start of parameter range
   * @param end Iterator pointing just beyond parameter range
   *
   * @return Reference to this
   */
  virtual Printer& beginFunc(const SizeVector::iterator& begin,
			     const SizeVector::iterator& end) = 0;
  /**
   * End of a function
   *
   * @return Reference to this
   */
  virtual Printer& endFunc() = 0;
  /**
   * Beginning of a term
   *
   * @param sign Sign character
   * @param fac  Factor
   *
   * @return Reference to this
   */
  virtual Printer& beginTerm(char sign, double fac)
  {
    mOut << "\n  " << sign << ' ';
    if (fac != 1) {
      mOut << fac;
      mult();
    }
    return *this;
  }
  /**
   * End of a term
   *
   * @return Reference to this
   */
  virtual Printer& endTerm() { mOut << std::flush; return *this; }
  /**
   * Print arguments
   *
   * @param begin Iterator pointing to start of argument range
   * @param end Iterator pointing just beyond argument range
   * @param sum If true, then print as a sum
   *
   * @return  Reference to this
   */
  virtual Printer& args(const SizeVector::iterator& begin,
			const SizeVector::iterator& end,
			bool sum=false)
  {
    char sep = sum ? '+' : ',';
    for (SizeVector::iterator i = begin; i != end; ++i) {
      if (i != begin) mOut << sep;
      mOut << 'h' << *i;
    }
    return *this;
  }
  /**
   * Print recursive call
   *
   * @param begin Iterator pointing to start of argument range
   * @param end Iterator pointing just beyond argument range
   *
   * @return  Reference to this
   */
  virtual Printer& call(const SizeVector::iterator& begin,
			const SizeVector::iterator& end)
  {
    Size n = end - begin;
    mOut << "QC" << n << '(';
    args(begin, end, false);
    mOut << ')' << std::flush;
    return *this;
  }
  /**
   * Print Q-vector access
   *
   * @param begin Iterator pointing to start of argument range
   * @param end Iterator pointing just beyond argument range
   *
   * @return  Reference to this
   */
  virtual Printer& access(const SizeVector::iterator& begin,
			  const SizeVector::iterator& end)
  {
    Size s = end - begin;
    // std::cerr << "Access with " << s << " elements to sum" << std::endl;
    mOut << "Q(";
    args(begin, end, true);
    mOut << ',' << s << ')' << std::flush;
    return *this;
  }
  /**
   * Print multiplication
   *
   * @return  Reference to this
   */
  virtual Printer& mult() { mOut << " * "; return *this; }
  /** output stream */
  std::ostream& mOut;
};

//____________________________________________________________________
struct Maxima : public Printer
{
  Maxima(std::ostream& out) : Printer(out) {}
  Printer& beginFunc(const SizeVector::iterator& begin,
			     const SizeVector::iterator& end)
  {
    Size n = end - begin;
    mOut << "QC" << n << '(';
    for (SizeVector::iterator i = begin; i != end; i++) {
      if (i != begin) mOut << ',';
      mOut << 'h' << *i;
    }
    mOut << ") :=";
    return *this;
  }
  Printer& endFunc() { mOut << ';' << std::endl; return *this; }
};

//____________________________________________________________________
struct Cxx : public Printer
{
  Cxx(std::ostream& out) : Printer(out) {}
  Printer& start()
  {
    mOut << "#include <complex>\n\n"
	 << "extern const std::complex<double>& Q(short h, unsigned short p);\n"
	 << std::endl;
    return *this;
  }
  Printer& stop(Size n)
  {
    SizeVector v(n);
    mOut << "std::complex<double> QC(const std::vector<short>& h)\n"
	 << "{\n"
	 << "  switch(h.size()) {\n";
    for (Size i = 1; i <= n; i++) {
      v[i-1] = i;
      mOut << "  case " << i << " return QC" << i << "(";
      args(v.begin(), v.begin()+i, false);
      mOut << ");\n";
    }
    mOut << "  }\n"
	 << "  return std::complex<double>(0,0);\n"
	 << "}\n" << std::endl;
    return *this;
  }

  Printer& beginFunc(const SizeVector::iterator& begin,
			     const SizeVector::iterator& end)
  {
    Size n = end - begin;
    mOut << "std::complex<double> QC" << n << '(';
    for (SizeVector::iterator i = begin; i != end; i++) {
      if (i != begin)
	mOut << ",\n                         ";
      mOut << "short h" << *i;
    }
    mOut << ")\n{\n  return";
    return *this;
  }
  Printer& endFunc() { mOut << ";\n}" << std::endl; return *this; }
  Printer& beginTerm(char sign, double fac)
  {
    mOut << "\n    " << sign << ' ';
    if (fac != 1) {
      mOut << fac << '.';
      mult();
    }
    return *this;
  }
};

//____________________________________________________________________
struct Latex : public Printer
{
  Latex(std::ostream& out) : Printer(out), mFirst(true) {}
  Printer& start()
  {
    mOut << "\\begin{align}\n";
    return *this;
  }
  Printer& stop(Size)
  {
    mOut << "\\end{align}" << std::endl;
    return *this;
  }
  Printer& beginFunc(const SizeVector::iterator& begin,
			     const SizeVector::iterator& end)
  {
    Size n = end - begin;
    mOut << "QC_{";
    args(begin, end, false);
    mOut << "}\\{" << n << "\\} &= ";
    mFirst = true;
    return *this;
  }
  Printer& endFunc()
  {
    mOut << "\\\\" << std::endl;
    return *this;
  }
  /**
   * Beginning of a term
   *
   * @param sign Sign character
   * @param fac  Factor
   *
   * @return Reference to this
   */
  Printer& beginTerm(char sign, double fac)
  {
    bool nl = !mFirst;
    mFirst  = false;
    if (nl) mOut << "\\nonumber\\\\";
    mOut << "\n  ";
    if (nl) mOut << " & ";
    mOut << sign << ' ';
    if (fac != 1) {
      mOut << fac;
      mult();
    }
    return *this;
  }
  /**
   * End of a term
   *
   * @return Reference to this
   */
  Printer& endTerm()
  {
    // mOut << "\\\\" << std::endl << "  & ";
    return *this;
  }
  Printer& call(const SizeVector::iterator& begin,
				const SizeVector::iterator& end)
  {
    Size n = end - begin;
    mOut << "QC_{";
    args(begin, end, false);
    mOut << "}\\{" << n << "\\}" << std::flush;
    return *this;
  }
  Printer& access(const SizeVector::iterator& begin,
			  const SizeVector::iterator& end)
  {
    Size s = end - begin;
    // std::cerr << "Access with " << s << " elements to sum" << std::endl;
    mOut << "Q_{";
    args(begin, end, true);
    mOut << "}^{" << s << '}' << std::flush;
    return *this;
  }
  Printer& mult() { mOut << " "; return *this; }
  bool mFirst;
};


//____________________________________________________________________
/**
 * Show the multi particle correlator as a recursive expression
 *
 * @param o    Output printer
 * @param n    Order of correlator
 * @param w    Current indices
 */
void _show_multi(Printer&           o,
		 const Size         n,
		 const SizeVector&  w,
		 SizeVector&        v)
{
  // std::cerr << "DoShowMulti(o,"<<n<<",w,\""<<fac<<"\");" << std::endl;
  if (n <= 1) {
    SizeVector a(1);
    a[0] = w[0];
    o.beginTerm('+',1).access(a.begin(), a.end()).endTerm(); return;
  }

  // Make vector of indices
  // SizeVector v(n-1);
  for (size_t i = 0; i < v.size(); i++) v[i] = w[i];

  Real    f = 1;
  Real    s = 1;
  Power   p = 1;
  for (Size m = n; m > 0; m--) {
    // Reset indices for combinations
    Size k  =  m-1;
    do {
      Real asf = std::abs(s*f);
      char sgn = (s * f >= 0 ? '+' : '-');

      o.beginTerm(sgn, asf);
      if (k > 0) o.call(v.begin(), v.begin()+k).mult();

      SizeVector u(v.begin()+k, v.begin()+n-1); u.push_back(w[n-1]);
      o.access(u.begin(), u.end());
    } while (next_combination(v.begin(), v.begin()+k, v.begin()+n-1));
    // Reset for next iteration
    for (size_t i = 0; i < v.size(); i++) v[i] = w[i];
    f       *= (n-k);
    s       *= -1;
    p++;
  }
}

//____________________________________________________________________
/**
 * Show the multi particle correlator as a recursive expression
 *
 * @param o    Output stream
 * @param n    Order of correlator
 * @param asMx As Maxima code
 */
void show_multi(const Size n,
		Printer& o)
{
  SizeVector v(n);
  SizeVector w(n);
  for (size_t i = 0; i < v.size(); i++) v[i] = i+1;

  o.beginFunc(v.begin(), v.end());

  _show_multi(o, n, v, w);
  o.endFunc();
}

//____________________________________________________________________
struct to_lower
{
  char operator()(char c) { return std::tolower(c); }
};

//____________________________________________________________________
void usage(const char* prog)
{
  std::cout << "Usage: " << prog << " [OPTIONS]\n\n"
	    << "Options:\n"
	    << "  -n  NUMBER   Largest correlator to print\n"
	    << "  -t  TYPE     How to print\n\n"
	    << "  -d           Show differential equations"
	    << "TYPE can be one of\n\n"
	    << "  maxima    Input to Maxima CAS\n"
	    << "  c++       C++ functions\n"
	    << "  latex     LaTeX equations\n"
	    << std::flush;
}

//____________________________________________________________________
int
main(int argc, char** argv)
{
  Size         n = 4;
  std::string  t = "maxima";
  // bool         d = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      std::cerr << argv[0] << ": Unknown option " << argv[i] << std::endl;
      return 1;
    }

    switch (argv[i][1]) {
    case 'h': usage(argv[0]); return 0;
    case 'n': n = atoi(argv[++i]); break;
    case 't': t = argv[++i]; break;
    // case 'd': d = true; break;
    default:
      std::cerr << argv[0] << ": Unknown option " << argv[i] << std::endl;
      return 1;
    }
  }

  std::transform(t.begin(), t.end(), t.begin(), to_lower());

  Printer* p = 0;
  if      (t == "maxima")            p = new Maxima(std::cout);
  else if (t == "cxx" || t == "c++") p = new Cxx(std::cout);
  else if (t == "latex")             p = new Latex(std::cout);
  if (!p) {
    std::cerr << "No printer defined for " << t << std::endl;
    return 1;
  }

  p->start();
  for (Size i = 1; i <= n; i++)
    show_multi(i, *p);
  p->stop(n);

  return 0;
}

//
// EOF
//
