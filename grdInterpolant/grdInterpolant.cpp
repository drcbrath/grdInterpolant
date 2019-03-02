// grdInterpolant.cpp

// version information
#define PROGNAME grdInterpolant
#define MAJOR 1
#define MINOR 0
#define PATCH 0
#define REV " "

#define STRX(s) #s
#define STR(s) STRX(x)

const int VerMahor = MAJOR;
const int VerMinor = MINOR;
const int VerPatch = PATCH;
const char VerRevision[] = REV;
const char VerString[] = STR(PROGNAME) " v" STR(MAJOR) "." STR(MINOR) "." STR(PATCH);

#include <iostream>
#include "grdInterpolant.h"

using namespace std;

#define CHECK_BIT(var,pos) ((var) & (li64<<(pos)))

int bplookup(std::vector<double> bp, double &x, int xtrpLt, int xtrpRt, size_t &nd)
{

}

grdInterpolant1D::grdInterpolant1D()
{
}

grdInterpolant1D::grdInterpolant1D(std::vector<double> xi, std::vector<double> yi)
{
}

grdInterpolant1D::~grdInterpolant1D()
{
}

double grdInterpolant1D::interp(double x1)
{
   size_t ndx;
   double y, ya, yb;
   double xa, xb;
   double s;

   // find index to breakpoint just below x
   if (bplookup(bp1, x1, xtrpLt, xtrpRt, ndx))
      return(std::numeric_limits<double>::quiet_NaN());

   // get breakpoints below and above x
   xa = bp1[ndx];
   xb = bp1[ndx + 1];
   ya = table[ndx];
   yb = table[ndx + 1];

   // linearly interpolate from table
   s = (x1 - xa) / (xb - xa);
   y = (1 - s)*ya + s * yb;

   return(y);
}


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

//------------------------------------------------------------------------------


//int lookup(double x)
//{
//	// bisection search for index of breakpoint element just below x,
//	//i.e. ndx such that Xs[ndx] <= x < Xs[ndx+1]
//	int na, nb, nc;
//
//	na = 0;
//	nb = this->size()-1;
//	nc = (nb - na) / 2;   // 1st guess in the middle
//
//	if(this->operator[](na) <= x && x < this->operator[](nb))
//	{
//		while (nb - na > 1)
//		{
//			if (this->operator[](nc) <= x)
//			{   // next guess in upper interval
//				na = nc;
//				nc = (nb - na) / 2;   
//			}
//			else
//			{   // next guess in lower interval
//				nb = nc;
//				nc = (nb - na) / 2;
//			}
//		}
//		return(na);
//	}
//	else if (this->operator[](nb) == x)
//	{
//		return(nb);
//	}
//
//	return(-1);
//
//}
