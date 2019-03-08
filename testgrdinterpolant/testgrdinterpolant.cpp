// testgrdinterpolant.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <math.h>
#include "grdInterpolant.h"

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

   // create breakpoints and table
   int N = 9;
   double pi = 3.14159236;
   double a;
   std::vector<double> bpX(N);
   std::vector<double> tableF(N);
   std::vector<double> tableP(N);
   for (int m = 0; m < N; m++) {
      // a = pi * ((n - 1) : -1 : 0)'/(n-1);
      a = pi * (N-1 - m) / (N-1);
      bpX[m] = cos(a);
      tableF[m] = sin(pi*bpX[m]);                             // nonlinear test function
      tableP[m] = cos(pi / 4)*bpX[m] + sin(pi / 4);   // linear test, therefore interpolation should be exact
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

   // create 1-D interpoplants
   grdInterpolant1D grd1DF(bpX, tableF);
   grdInterpolant1D grd1DP(bpX, tableP);

   // test
   std::cout << "\n1-D -- extrapolation LINEAR\n";
   std::cout << "x << p << pe << p - pe << f << fe << f - fe" << std::endl;
   double f, p, fe, pe;
   for (int m = 0; m < x.size(); m++) {
      f = grd1DF.interp(x[m]);
      p = grd1DP.interp(x[m]);
      fe = sin(pi*x[m]);
      pe = cos(pi / 4)*x[m] + sin(pi / 4);
      // print comparison
      std::cout << x[m] << " " << p << " " << pe << " " << p - pe << " " << f << " " << fe << " " << f - fe << std::endl;
   }

   grd1DF.xtrpLt = NEAREST; grd1DF.xtrpRt = NEAREST;
   grd1DP.xtrpLt = NEAREST; grd1DP.xtrpRt = NEAREST;
   std::cout << "\n1-D -- extrapolation NEAREST, left & right\n";
   std::cout << "x << p << pe << p - pe << f << fe << f - fe" << std::endl;
   for (int m = 0; m < x.size(); m++) {
      f = grd1DF.interp(x[m]);
      p = grd1DP.interp(x[m]);
      fe = sin(pi*x[m]);
      pe = cos(pi / 4)*x[m] + sin(pi / 4);
      // print comparison
      std::cout << x[m] << " " << p << " " << pe << " " << p - pe << " " << f << " " << fe << " " << f - fe << std::endl;
   }

   grd1DF.xtrpLt = NONE; grd1DF.xtrpRt = NONE;
   grd1DP.xtrpLt = NONE; grd1DP.xtrpRt = NONE;
   std::cout << "\n1-D -- extrapolation NONE, left & right\n";
   std::cout << "x << p << pe << p - pe << f << fe << f - fe" << std::endl;
   for (int m = 0; m < x.size(); m++) {
      f = grd1DF.interp(x[m]);
      p = grd1DP.interp(x[m]);
      fe = sin(pi*x[m]);
      pe = cos(pi / 4)*x[m] + sin(pi / 4);
      // print comparison
      std::cout << x[m] << " " << p << " " << pe << " " << p - pe << " " << f << " " << fe << " " << f - fe << std::endl;
   }

   // 2-D

   // 3-D

   // N-D

}
