#include <iostream>
#include <fstream>

#include "bsplines.h"

typedef unsigned int uint;

uint pingcount = 0;

void ping()
{
	cout << "ping" << pingcount << endl;
	pingcount++;
}

int main ()
{
	string knotfile = "testknots.txt";
	string gridfile = "testgrid.txt";

	int order = 4;
	double tol = 1e-9;

	int xmin = -2;
	int xmax = 2;
	int ksteps = 5;
	double gsteps = 100;

	ofstream kstream, gstream;
	stream.open(knotfile);
	gstream.open(gridfile);

	for (int k = 0; k < ksteps; k++)
	{
		kstream << xmin + k << endl;
	}

	double dg = (xmax - xmin) / gsteps;
	for (int k = 0; k <= gsteps; k++)
	{
		double x = xmin + dg*k;

		gstream << x << endl;
	}

	kstream.close();
	gstream.close();

	bsplines splines (knotfile,"",order,tol);

	splines.calcspline(-1);
	//splines.calcsplines();

	splines.writesplines();

	ping();
}
