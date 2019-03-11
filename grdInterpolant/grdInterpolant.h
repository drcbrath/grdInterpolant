#pragma once

#include <array>
#include <vector>

//------------------------------------------------------------------------------

enum xtraps { NONE, NEAREST, LINEAR };

class grdInterpolant1D
{
public:
   // extrapolation flags; default linear on left and right end of breakpoint vectors
   int xtrpLt = LINEAR;
   int xtrpRt = LINEAR;

   std::vector<double> bp1;      // vector of breakpoints for 1st independent variable
   std::vector<double> table;    // 1-D table as a vector, values to interpolate

   grdInterpolant1D();
   grdInterpolant1D(std::vector<double> xi, std::vector<double> yi);

   ~grdInterpolant1D();

   double interp(double x);

private:

};

class grdInterpolant2D
{
public:
   // extrapolation flags; default linear on left and right end of breakpoint vectors
   std::vector<int> xtrpLt = { LINEAR, LINEAR };
   std::vector<int> xtrpRt = { LINEAR, LINEAR };

   std::vector<double> bp1;                     // vector of breakpoints for independent variable 1
   std::vector<double> bp2;                     // vector of breakpoints for independent variable 2
   std::vector< std::vector<double> > table;    // 2-D table as a vector of vectors, values to interpolate

   grdInterpolant2D();
   grdInterpolant2D(std::vector<double> xi, std::vector<double> yj, std::vector< std::vector<double> > vij);

   ~grdInterpolant2D();

   double interp(double x1, double x2);

private:

};


class grdInterpolant3D
{
public:
   // extrapolation flags; default linear on left and right end of breakpoint vectors
   std::vector<int> xtrpLt = { LINEAR, LINEAR, LINEAR };
   std::vector<int> xtrpRt = { LINEAR, LINEAR, LINEAR };

   std::vector<double> bp1;      // vector of breakpoints for independent variable 1
   std::vector<double> bp2;      // vector of breakpoints for independent variable 2
   std::vector<double> bp3;      // vector of breakpoints for independent variable 4
   std::vector< std::vector< std::vector<double> > > table;    // 3-D table as a vector, values to interpolate

   grdInterpolant3D();
   grdInterpolant3D(std::vector<double>xi, std::vector<double>xj, std::vector<double>yk, std::vector< std::vector< std::vector<double> > > vijk);

   ~grdInterpolant3D();

   double interp(double x1, double x2, double x3);

private:

};


class grdInterpolantND
{
public:
   size_t nDims;   // number of dimensions, aka tensor rank
   size_t nverts;   // number of table vertices surrounding inputpoint (x1, x2, x3, ..., xn) at which to interpolate
   std::vector<size_t> sizes; // vector of size of each dimension (like return from matlab size(table) )
   std::vector<int> xtrpLt;
   std::vector<int> xtrpRt;

   std::vector< std::vector<double> > bps;   // vector of breakpoint vectors for independent variables
   std::vector<double> table;   // N-D table as a vector (needs indexing function to access like a multi-dimensional table)

   grdInterpolantND();
   grdInterpolantND(std::vector< std::vector<double> > bpXs, std::vector<double> bpValues);

   ~grdInterpolantND();

   double interp(std::vector<double> x);

private:
   std::vector<size_t> Nofs;   // vector of offsets from subscript indices into table(j1,j2,j3, ..., jn) to linear index table(i)
   std::vector<size_t> nds;    // vector of indices found from lookup of input x into each breakpoint vector
   std::vector<double> xf;     // vector of graction of distance of xi between dim[m] endpoints; i.e. xa --- xi ------xb, xf = dist(xa,xi)/dist(xa,xb)
   std::vector<double> cxf;    // vector of 1-xf

};
