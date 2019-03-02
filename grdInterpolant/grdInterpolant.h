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

grdInterpolant1D::grdInterpolant1D()
{
}

grdInterpolant1D::grdInterpolant1D(std::vector<double> xi, std::vector<double> yi)
{
}

grdInterpolant1D::~grdInterpolant1D()
{
}

//double grdInterpolant1D::interp(double x)
//{
//}

class grdInterpolant2D
{
public:
   // extrapolation flags; default linear on left and right end of breakpoint vectors
   std::vector<int> xtrpLt;
   std::vector<int> xtrpRt;

   std::vector<double> bp1;      // vector of breakpoints for independent variable 1
   std::vector<double> bp2;      // vector of breakpoints for independent variable 2
   std::vector<double> table;    // 2-D table as a vector, values to interpolate

   grdInterpolant2D();
   grdInterpolant2D(std::vector<double> bp1, std::vector<double> bp2, std::vector< std::vector<double> > table);

   ~grdInterpolant2D();

   double interp(double x1, double x2);

private:

};

grdInterpolant2D::grdInterpolant2D()
{
}

grdInterpolant2D::grdInterpolant2D(std::vector<double> bp1, std::vector<double> bp2, std::vector< std::vector<double> > table)
{
}

grdInterpolant2D::~grdInterpolant2D()
{
}

double grdInterpolant2D::interp(double x1, double x2)
{
   return 0.0;
}

class grdInterpolant3D
{
public:
   // extrapolation flags; default linear on left and right end of breakpoint vectors
   std::vector<int> xtrpLt;
   std::vector<int> xtrpRt;

   std::vector<double> bp1;      // vector of breakpoints for independent variable 1
   std::vector<double> bp2;      // vector of breakpoints for independent variable 2
   std::vector<double> bp3;      // vector of breakpoints for independent variable 4
   std::vector<double> table;    // 3-D table as a vector, values to interpolate

   grdInterpolant3D();
   grdInterpolant3D(std::vector<double>bp1, std::vector<double>bp2, std::vector<double>bp3, std::vector< std::vector< std::vector<double> > >table);

   ~grdInterpolant3D();

   double interp(double x1, double x2, double x3);

private:

};

grdInterpolant3D::grdInterpolant3D()
{
}

grdInterpolant3D::grdInterpolant3D(std::vector<double>bp1, std::vector<double>bp2, std::vector<double>bp3, std::vector< std::vector< std::vector<double> > >table)
{
}

grdInterpolant3D::~grdInterpolant3D()
{
}

double grdInterpolant3D::interp(double x1, double x2, double x3)
{
   return 0.0;
}

class grdInterpolantND
{
public:
   size_t nDims;   // number of dimensions, aka tensor rank
   size_t nverts;   // number of table vertices surrounding inputpoint (x1, x2, x3, ..., xn) at which to interpolate
   std::vector<size_t> sizes; // vector of size of each dimension (like return from matlab size(table) )
   std::vector<int> xtrpLt;
   std::vector<int> xtrpRt;

   std::vector< std::array<char, 2> > exfs;   // vector of extrapolation flags ?
   std::vector< std::vector<double> > bps;   // vector of breakpoint vectors for independent variables

   std::vector<double> table;   // N-D table as a vector (needs indexing function to access like a multi-dimensional table)

   grdInterpolantND();
   grdInterpolantND(std::vector< std::vector<double> > bps, std::vector<double> table);

   ~grdInterpolantND();

   double interp(std::vector<double> x);    // handles singleton dimensions
   double interps(std::vector<double> x);   // for speed, but does not handle singleton dimensions

private:
   std::vector<size_t> Nofs;   // vector of offsets from subscript indices into table(j1,j2,j3, ..., jn) to linear index table(i)
   std::vector<size_t> nds;    // vector of indices found from lookup of input x into each breakpoint vector
   std::vector<double> xf;     // vector of graction of distance of xi between dim[m] endpoints; i.e. xa --- xi ------xb, xf = dist(xa,xi)/dist(xa,xb)
   std::vector<double> cxf;    // vector of 1-xf

};

grdInterpolantND::grdInterpolantND()
{
}

grdInterpolantND::grdInterpolantND(std::vector<std::vector<double>> bps, std::vector<double> table)
{
}

grdInterpolantND::~grdInterpolantND()
{
}

double grdInterpolantND::interp(std::vector<double> x)
{
   return 0.0;
}

double grdInterpolantND::interps(std::vector<double> x)
{
   return 0.0;
}
