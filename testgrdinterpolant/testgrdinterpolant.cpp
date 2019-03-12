// testgrdInterpolant -- test grdInterpolant class, interpolation of gridded data
//
// Copyright(c) 2012 - 2019 Curtis B.Rath, PhD., aka drcbrath
// available under the MIT license from Github, https://github.com/drcbrath/mdlAtmos

#include "pch.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include "grdInterpolant.h"

inline double F2(double x, double y)
{
   return(cos(M_PI*x)*sin(M_PI*y));
}

inline double F3(double x, double y, double z)
{
   return(cos(M_PI*x)*sin(M_PI*y) + tan(M_PI*z));
}

int main()
{
   // Test grdInterpolant class
   //
   // 1-D
   //
   // lookup, extrapolation handling, and 1-D interpolation
   // lookup & extrapolation handling are identical (algorithm & code) for all
   // dimensions (1-D, 2-D, N-D), so it suffices to test once here
   //

   // need to make sure bpX is strictly increasing, else lookup almost always wrong

   std::cout << "Test grdInterpolant\n";

   // create breakpoints and tables
   int N = 9;
   double a;
   std::vector<double> bpX(N);
   std::vector<double> bpY(N);
   std::vector<double> bpZ(N);
   std::vector<double> tableF1(N);
   std::vector<double> tableP1(N);
   std::vector<std::vector<double> > tableF2(N, std::vector<double>(N));
   std::vector<std::vector<double> > tableP2(N, std::vector<double>(N));
   std::vector< std::vector<std::vector<double> > > tableF3(N, tableF2);
   std::vector< std::vector<std::vector<double> > > tableP3(N, tableP2);

   double df1dx, dfdx, dfdy, dfdz, g1, p1, p2, p3, f10, f20, f30;
   double x0 = 0.25, y0 = 0.25, z0 = 0.25;

   double xt = 0.25, yt = 0.25, zt = 0.25;

   f10 = cos(M_PI*x0);
   f20 = F2(x0,y0);
   f30 = F3(x0,y0,z0);
   df1dx = -sin(M_PI*x0);
   dfdx = -M_PI * sin(M_PI*x0)*sin(M_PI*y0);
   dfdy = M_PI * cos(M_PI*x0)*cos(M_PI*y0);
   dfdz = M_PI * 1.0/pow(cos(M_PI*z0),2);
   p1 = df1dx * (xt - x0) + f10;
   p2 = dfdx * (xt - x0) + dfdy * (yt - y0) + f20;
   p3 = dfdx * (xt - x0) + dfdy * (yt - y0) + dfdz * (zt - z0) + f30;

   for (int k = 0; k < N; k++) {
      a = M_PI * (N - 1 - k) / (N - 1);
      bpX[k] = cos(a);
      bpY[k] = cos(a);
      bpZ[k] = cos(a);
   }
    
   for (int k = 0; k < N; k++) {
      tableF1[k] = sin(M_PI*bpX[k]);                     // nonlinear test function
      tableP1[k] = df1dx*(bpX[k]-x0) + f10;            // linear test, therefore interpolation should be exact
      for (int j = 0; j < N; j++)
      {
         tableF2[k][j] = F2(bpX[j],bpY[k]);
         tableP2[k][j] = dfdx * (bpX[j] - x0) + dfdy * (bpY[k] - y0) + f20;
         for (int i = 0; i < N; i++)
         {
            tableF3[k][j][i] = F3(bpX[i], bpY[j], bpZ[k]);
            tableP3[k][j][i] = dfdx * (bpX[i] - x0) + dfdy * (bpY[j] - y0) + dfdz * (bpZ[k] - z0) + f30;
         }
      }
   }

   // test interpolation points
   std::vector<double> x = {bpX[0] - 0.5,                                    // below table
                            bpX[0],                                          // on 1st breakpoint
                            0.5*(bpX[0] + bpX[1]),                  // in 1st interval
                            bpX[1],                                          // right end of 1st interval
                            bpX[N / 2],                                      // on a middle breakpoint
                            0.5*(bpX[N / 2] + bpX[N / 2 + 1]),  // in a middle interval   
                            bpX[N - 2],                                      // on next to last breakpoint
                            0.5*(bpX[N - 2] + bpX[N - 1]),      // in last interval
                            bpX[N - 1],                                      // on last point
                            bpX[N - 1] + 0.5};                               // above table

   // create interpolants
   grdInterpolant1D grd1DF(bpX, tableF1);
   grdInterpolant1D grd1DP(bpX, tableP1);

   // test
   std::cout << "\n1-D -- extrapolation LINEAR\n";
   std::cout << "x << p << pe << p - pe << f << fe << f - fe" << std::endl;
   double p, fe, pe;
   for (int m = 0; m < x.size(); m++) {
      g1 = grd1DF.interp(x[m]);
      p = grd1DP.interp(x[m]);
      fe = sin(M_PI*x[m]);
      pe = cos(M_PI / 4)*x[m] + sin(M_PI / 4);
      // print comparison
      std::cout << x[m] << " " << p << " " << pe << " " << p - pe << " " << g1 << " " << fe << " " << g1 - fe << std::endl;
   }

   grd1DF.xtrpLt = NEAREST; grd1DF.xtrpRt = NEAREST;
   grd1DP.xtrpLt = NEAREST; grd1DP.xtrpRt = NEAREST;
   std::cout << "\n1-D -- extrapolation NEAREST, left & right\n";
   std::cout << "x << p << pe << p - pe << f << fe << f - fe" << std::endl;
   for (int m = 0; m < x.size(); m++) {
      g1 = grd1DF.interp(x[m]);
      p = grd1DP.interp(x[m]);
      fe = sin(M_PI*x[m]);
      pe = cos(M_PI / 4)*x[m] + sin(M_PI / 4);
      // print comparison
      std::cout << x[m] << " " << p << " " << pe << " " << p - pe << " " << g1 << " " << fe << " " << g1 - fe << std::endl;
   }

   grd1DF.xtrpLt = NONE; grd1DF.xtrpRt = NONE;
   grd1DP.xtrpLt = NONE; grd1DP.xtrpRt = NONE;
   std::cout << "\n1-D -- extrapolation NONE, left & right\n";
   std::cout << "x << p << pe << p - pe << f << fe << f - fe" << std::endl;
   for (int m = 0; m < x.size(); m++) {
      g1 = grd1DF.interp(x[m]);
      p = grd1DP.interp(x[m]);
      fe = sin(M_PI*x[m]);
      pe = cos(M_PI / 4)*x[m] + sin(M_PI / 4);
      // print comparison
      std::cout << x[m] << " " << p << " " << pe << " " << p - pe << " " << g1 << " " << fe << " " << g1 - fe << std::endl;
   }

   // 2-D
   // create interpolants
   grdInterpolant2D grd2DF(bpX, bpY, tableF2);
   grdInterpolant2D grd2DP(bpX, bpY, tableP2);

   double f2, gf2, gp2;

   f2 = F2(xt, yt);              // exact
   gf2 = grd2DF.interp(xt, yt);  // interpolated
   gp2 = grd2DP.interp(xt, yt);

   std::cout << std::endl << xt << " " << yt << " " << gp2 << " " << p2 << " " << gp2 - p2 << " " << gf2 << " " << f2 << " " << gf2 - f2 << std::endl;

   // 3-D
   grdInterpolant3D grd3DF(bpX, bpY, bpY, tableF3);
   grdInterpolant3D grd3DP(bpX, bpY, bpZ, tableP3);

   double f3, gf3, gp3;

   f3 = F3(xt, yt, zt);              // exact
   gf3 = grd3DF.interp(xt, yt, zt);  // interpolated
   gp3 = grd3DP.interp(xt, yt, zt);

   std::cout << std::endl << xt << " " << yt << " " << zt << " " << gp3 << " " << p3 << " " << gp3 - p3 << " " << gf3 << " " << f3 << " " << gf3 - f3 << std::endl;

   // N-D
   // create vector of vectors for breakpoints
   // create table as vector
   // test with arbitrary dimension, say 6
   size_t nDims = 6;
   std::vector<std::vector<double>> bpXs(nDims);   // vector of vectors of breakpoints
   std::vector<double> dpdx(nDims);   // coefficients to construct linear test function
   std::vector<double> x00 = {0.25,0.25,0.25,0.25,0.25,0.25};   // anchor point of linear test function

   double f0=0;

   for (size_t j = 0; j < nDims; j++)
   {
      f0 += sin(M_PI*x00[j]);
   }
   for (size_t j = 0; j < nDims; j++)
   {
      dpdx[j] = M_PI * cos(M_PI*x00[j]) + f0;
   }
   double pn = f0;

   std::vector<size_t> insizes = { 9,9,9,9,9,9 };   // can change size of each dimensions breakpoints here, smaller for faster test; but worse FN approximation
   for (size_t j = 0; j < nDims; j++)
   {
      bpXs[j].resize(insizes[j]);
      for (size_t i = 0; i < bpXs[j].size(); i++)
      {
         a = M_PI * (N - 1 - i) / (N - 1);
         bpXs[j][i] = cos(a);
      }
   }

   // compute number of elements in the table
   size_t nTableNumel = 1;
   for (size_t i = 0; i < nDims; i++)
   {
      nTableNumel *= bpXs[i].size();
   }
   std::vector<double> tableFN(nTableNumel);   // create table with required number of elements
   std::vector<double> tablePN(nTableNumel);   // create table with required number of elements

   // fill the table
   int n = 0;   // initialize table element index
   std::vector<int> nds(6);   // breakpoint vector indices for each dimension
   for (nds[5] = 0; nds[5] < bpXs[5].size(); nds[5]++) {
      for (nds[4] = 0; nds[4] < bpXs[4].size(); nds[4]++) {
         for (nds[3] = 0; nds[3] < bpXs[3].size(); nds[3]++) {
            for (nds[2] = 0; nds[2] < bpXs[2].size(); nds[2]++) {
               for (nds[1] = 0; nds[1] < bpXs[1].size(); nds[1]++) {
                  for (nds[0]=0; nds[0] < bpXs[0].size(); nds[0]++){
                     double vf = 0, vp = f0;
                     for (size_t j = 0; j < nDims; j++)
                     {
                        vf += sin(M_PI*bpXs[j][nds[j]]);
                        vp += M_PI * dpdx[j] * (bpXs[j][nds[j]] - x00[j]);
                     }
                     tableFN[n] = vf;
                     tablePN[n] = vp;
                     n++;   // increment table element index
                  }
               }
            }
         }
      }
   }
   
   grdInterpolantND grdFND(bpXs, tableFN);
   grdInterpolantND grdPND(bpXs, tablePN);

   double gfn, gpn;
   gfn = grdFND.interp(x00);
   gpn = grdPND.interp(x00);

   std::cout << std::endl << "N-D test" << std::endl;
   std::cout << std::endl << gpn << " " << pn << " " << gpn - pn << " " << gfn << " " << f0 << " " << gfn - f0 << std::endl;

}
