// grdInterpolant.cpp

// version information
#define PROGNAME grdInterpolant
#define MAJOR 1
#define MINOR 0
#define PATCH 0
#define REV " "

#define STRX(s) #s
#define STR(s) STRX(x)

const int VerMajor = MAJOR;
const int VerMinor = MINOR;
const int VerPatch = PATCH;
const char VerRevision[] = REV;
const char VerString[] = STR(PROGNAME) " v" STR(MAJOR) "." STR(MINOR) "." STR(PATCH);

#include "grdInterpolant.h"

using namespace std;

#define CHECK_BIT(var,pos) ((var) & (1i64<<(pos)))

int bplookup(std::vector<double> bp, double &x, int xtrpLt, int xtrpRt, size_t &nd)
{

   	// bisection search for index of breakpoint element just below x,
   	//i.e. ndx such that Xs[ndx] <= x < Xs[ndx+1]
   	int na, nb, nc;

      if (isnan(x))
         return(-1);   // fail, cannot interpolate at NaN                                                   

      na = 0;
      nb = bp.size() - 1;

      if (x < bp[0])  // x below interval (note: okay for singleton dimension)
      {
         nd = 0;      // set to use 1st subinterval
         // handle per left extrapolation flag
         switch (xtrpLt)
         {
         case NONE:
            return(-1);   // fail
            break;
         case NEAREST:
            x = bp[0];   // bp nearest to x, so extrapolation will evaluate to value at nearest bp
            break;
         case LINEAR:
         default:
            break;
         }
         return(0);   // success
      }
      else if (bp[na] <= x && x < bp[nb])   // x within interval (note: will fail for singleton dimension, as desired)
      {
            while (nb - na > 1)
            {
               nc = (nb + na) / 2;   // guess in the middle
               if (bp[nc] <= x)      // next guess in upper interval
                  na = nc;
               else if (bp[nc] == x)
               {                     // found it
                  na = nc;
                  break;
               }
               else                  // next guess in lower interval
                  nb = nc;
            }
            nd = na;     // found it
            return(0);   // return success
         }
      else if (bp[nb] == x)      // found it on last breakpoint
      {
         nd = nb-1 ;
         return(0);   // return success
      }
      else   // x above interval (note: okay for singleton dimension)
      {
         nd = nb - 1;   // set to use last subinterval
         // handle per right extrapolation flag
         switch (xtrpRt)
         {
         case NONE:
            x = NAN;
            return(-1);   // fail
            break;
         case NEAREST:
            x = bp[nb];   // bp nearest to x, so extrapolation will evaluate to value at nearest bp
            break;
         case LINEAR:
         default:
            break;
         }
         return(0);   // success
      }
   
   	return(-1);   // Now, how did this happen?! return failure

}

grdInterpolant1D::grdInterpolant1D()
{
}

grdInterpolant1D::grdInterpolant1D(std::vector<double> xi, std::vector<double> yi)
{
   bp1 = xi;
   table = yi;
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

//------------------------------------------------------------------------------

grdInterpolant2D::grdInterpolant2D()
{
}

grdInterpolant2D::grdInterpolant2D(std::vector<double> xi, std::vector<double> yj, std::vector< std::vector<double> > vij)
{
   bp1 = xi;
   bp2 = yj;
   table = vij;
}

grdInterpolant2D::~grdInterpolant2D()
{
}

double grdInterpolant2D::interp(double x1, double x2)
{
   size_t nd1, nd2;
   double x1a, x1b, x2a, x2b;
   double v00, v01, v10, v11;
   double s, t, v;

   // find index to breakpoint just below x
   if (bplookup(bp1, x1, xtrpLt[0], xtrpRt[0], nd1))
      return(std::numeric_limits<double>::quiet_NaN());
   if (bplookup(bp2, x2, xtrpLt[1], xtrpRt[1], nd2))
      return(std::numeric_limits<double>::quiet_NaN());

   // get corner points of rectangle to interpolate within
   x1a = bp1[nd1];
   x1b = bp1[nd1 + 1];
   x2a = bp2[nd2];
   x2b = bp2[nd2 + 1];
   // get table values at corner points
   v00 = table[nd2][nd1];
   v01 = table[nd2][nd1 + 1];
   v10 = table[nd2+1][nd1];
   v11 = table[nd2+1][nd1 + 1];

   // x as fractions of distance from left (a) to right (b) ends in each dimension
   s = (x1 - x1a) / (x1b - x1a);
   t = (x2 - x2a) / (x2b - x2a);

   // interpolate
   v = (1 - t)*((1 - s)*v00 + s * v01) + t * ((1 - s)*v10 + s * v11);

   return(v);
}

//------------------------------------------------------------------------------

grdInterpolant3D::grdInterpolant3D()
{
}

grdInterpolant3D::grdInterpolant3D(std::vector<double>xi, std::vector<double>yj, std::vector<double>zk, std::vector< std::vector< std::vector<double> > > vijk)
{
   bp1 = xi;
   bp2 = yj;
   bp3 = zk;
   table = vijk;
}

grdInterpolant3D::~grdInterpolant3D()
{
}

double grdInterpolant3D::interp(double x1, double x2, double x3)
{
   size_t nd1, nd2, nd3;
   double x1a, x1b, x2a, x2b, x3a, x3b;
   double v000, v001, v010, v011, v100, v101, v110, v111;
   double s, t, u, v;

   // find index to breakpoint just below x
   if (bplookup(bp1, x1, xtrpLt[0], xtrpRt[0], nd1))
      return(std::numeric_limits<double>::quiet_NaN());
   if (bplookup(bp2, x2, xtrpLt[1], xtrpRt[1], nd2))
      return(std::numeric_limits<double>::quiet_NaN());
   if (bplookup(bp3, x3, xtrpLt[2], xtrpRt[2], nd3))
      return(std::numeric_limits<double>::quiet_NaN());

   // get corner points of cell to interpolate within
   x1a = bp1[nd1];
   x1b = bp1[nd1 + 1];
   x2a = bp2[nd2];
   x2b = bp2[nd2 + 1];
   x3a = bp3[nd3];
   x3b = bp3[nd3 + 1];
   // get table values at corner points
   v000 = table[nd3][nd2][nd1];
   v001 = table[nd3][nd2][nd1+1];
   v010 = table[nd3][nd2+1][nd1];
   v011 = table[nd3][nd2+1][nd1+1];
   v100 = table[nd3+1][nd2][nd1];
   v101 = table[nd3+1][nd2][nd1+1];
   v110 = table[nd3+1][nd2+1][nd1];
   v111 = table[nd3+1][nd2+1][nd1+1];

   // x as fractions of distance from left (a) to right (b) ends in each dimension
   s = (x1 - x1a) / (x1b - x1a);
   t = (x2 - x2a) / (x2b - x2a);
   u = (x3 - x3a) / (x3b - x3a);

   v = (1 - u)*(1 - t)*(1 - s)*v000 + (1 - u)*(1 - t)*(s)*v001 +
       (1 - u)*(t)*(1 - s)*v010 + (1 - u)*(t)*(s)*v011 +
       (u)*(1 - t)*(1 - s)*v100 + (u)*(1 - t)*(s)*v101 +
       (u)*(t)*(1 - s)*v110 + (u)*(t)*(s)*v111;

   return(v);
}

//------------------------------------------------------------------------------

grdInterpolantND::grdInterpolantND()
{
}

grdInterpolantND::grdInterpolantND(std::vector< std::vector<double> > bpXs, std::vector<double> bpValues)
{
   bps = bpXs;         // copy vector of breakpoint vectors
   table = bpValues;   // copy table values at breakpoints

   nDims = bps.size();       // number of dimensions, aka tensor rank
   nverts = 1i64 << nDims;   // pow(2, nDims);   // number of table vertices surrounding input point (x1, x2, x3, ..., xn) at which to interpolate

   // resize per nDims to make room
   sizes.resize(nDims);
   xtrpLt.resize(nDims);
   xtrpRt.resize(nDims);
   Nofs.resize(nDims);   // vector of offsets from subscript indices into table(j1,j2,j3, ..., jn) to linear index table(i)
   nds.resize(nDims);    // vector of indices found from lookup of input x into each breakpoint vector
   xf.resize(nDims);     // vector of graction of distance of xi between dim[m] endpoints; i.e. xa --- xi ------xb, xf = dist(xa,xi)/dist(xa,xb)
   cxf.resize(nDims);    // vector of 1-xf

  for (size_t i = 0; i < nDims; i++)   
   {
      sizes[i] = bps[i].size();    // vector of size of each dimension (like return from matlab size(table) )
      xtrpLt[i] = LINEAR;          // default extrapolation handling
      xtrpRt[i] = LINEAR;
   }

  Nofs[0] = 1;   // initialize cummulative product
  for (size_t i = 0; i < nDims-1; i++)
  {
     Nofs[i + 1] = sizes[i] * Nofs[i];
  }
   
}

grdInterpolantND::~grdInterpolantND()
{
}

double grdInterpolantND::interp(std::vector<double> x)
{
   size_t mask;
   size_t Nxiv;
   double v;

  size_t Nxi = 0;
  for (size_t m = 0; m < nDims; m++)
   {
      if (bplookup(bps[m], x[m], xtrpLt[m], xtrpRt[m], nds[m]))
         return(std::numeric_limits<double>::quiet_NaN());

      if (sizes[m] != 1)   // not a singleton dimension
      {
         xf[m] = (x[m] - bps[m][nds[m]]) / (bps[m][nds[m] + 1] - bps[m][nds[m]]);
         cxf[m] = 1 - xf[m];
      }
      else   // is a singleton dimension
      {
         xf[m] = 0.0;
         cxf[m] = 1.0;
      }

      Nxi += (nds[m])*Nofs[m];
   }

   double vi = 0.0;
   for (size_t n = 0; n < nverts; n++)
   {
      Nxiv = Nxi;
      for (size_t m = 0; m < nDims; m++)
      {
         if (CHECK_BIT(n, m)) Nxiv += Nofs[m];
      }
      v = table[Nxiv];
      for (size_t m = 0; m < nDims; m++)
      {
         mask = CHECK_BIT(n, m);
         if (mask == 0)
            v = v * cxf[m];
         else
            v = v * xf[m];
      }
      vi += v;
   }
   return(vi);
}
