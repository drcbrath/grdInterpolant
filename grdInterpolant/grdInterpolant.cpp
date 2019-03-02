// grdInterpolant.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <iostream>

using namespace std;

class bpVector: public vector<double>
{
public:
	bpVector();
	~bpVector();

	int lookup(double x);

private:

};

bpVector::bpVector()
{
}

bpVector::~bpVector()
{
}

int bpVector::lookup(double x)
{
	// bisection search for index of breakpoint element just below x,
	//i.e. ndx such that Xs[ndx] <= x < Xs[ndx+1]
	int na, nb, nc;

	na = 0;
	nb = this->size()-1;
	nc = (nb - na) / 2;   // 1st guess in the middle

	if(this->operator[](na) <= x && x < this->operator[](nb))
	{
		while (nb - na > 1)
		{
			if (this->operator[](nc) <= x)
			{   // next guess in upper interval
				na = nc;
				nc = (nb - na) / 2;   
			}
			else
			{   // next guess in lower interval
				nb = nc;
				nc = (nb - na) / 2;
			}
		}
		return(na);
	}
	else if (this->operator[](nb) == x)
	{
		return(nb);
	}

	return(-1);

}

class grdInterpolant1D
{
public:

	bpVector bp;
	std::vector<double> table;

	grdInterpolant1D();
	grdInterpolant1D(int nx);
	grdInterpolant1D(vector<double> xi, vector<double> yi);

	~grdInterpolant1D();

	double Interp(double x);

private:

};

grdInterpolant1D::grdInterpolant1D()
{
}

grdInterpolant1D::grdInterpolant1D(int nx)
{
	this->bp.resize(nx);
	this->table.resize(nx);
}

grdInterpolant1D::grdInterpolant1D(vector<double> xi, vector<double> yi)
{
	bp.resize(xi.size());
	for (int i = 0; i < bp.size(); i++)
	{
		bp[i] = xi[i];
	}
	table = yi;
}

grdInterpolant1D::~grdInterpolant1D()
{
}

double grdInterpolant1D::Interp(double x)
{
	int ndx;
	double y, ya, yb;
	double xa, xb;
	double s;

	// find index to breakpoint just below x
	ndx = bp.lookup(x);

	// get breakpoints below and above x
	xa = bp[ndx];
	xb = bp[ndx + 1];
	ya = table[ndx];
	yb = table[ndx + 1];

	// linearly interpolate from table
	s = (x - xa) / (xb - xa);
	y = (1 - s)*ya + s*yb;

	return(y);
}

int main()
{
	int nx = 5;
	double xi[] = { 0.0, 1.1, 2.2, 3.3, 4.4 };
	double yi[] = { 0.0, 1.21, 4.84, 10.89, 17.76 };

	vector<double> xv = { 0.0, 1.1, 2.2, 3.3, 4.4 };
	vector<double> yv = { 0.0, 1.21, 4.84, 10.89, 17.76 };

	double x;

	grdInterpolant1D gi(5);

	for (int i = 0; i <= nx - 1; i++)
	{
		gi.bp[i] = xi[i];
		gi.table[i] = yi[i];
	}

	x = gi.Interp(0.55);
	cout << x << endl;

	grdInterpolant1D gv(xv, yv);

	return 0;
}