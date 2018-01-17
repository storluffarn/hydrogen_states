
#include <cmath>
#include <iostream>
#include <fstream>
#include <armadillo>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include "bsplines.h"

using namespace std;
typedef unsigned int uint;

double pi = 4*atan(1);

int pingcount = 0;
void ping() {cout << "ping" << pingcount << endl; pingcount++;}

constexpr uint64_t factorial(uint64_t n)
{ 
	return n == 0 ? 1 : n*factorial(n-1); 
}

double nlegendre (uint n)
{
	double norm = sqrt(2/(2*n+1));

	return norm;
}

double dlegendre (double x, uint n)
{
	double accdled = 0;

	for (int k = n; k >= 0; k--)
	{
		double t = 2*boost::math::legendre_p(k,x)/(nlegendre(k));

		accdled += t;
	}

	return accdled;
}

double weight (double x, uint n)
{
	double w = 2/((1-pow(x,2))*pow(dlegendre(x,n),2));
	
	return w;
}

void makeknots (double rmin, double rmax, uint steps, uint mode, string filename)
{
	vector <double> grid;
	grid.reserve(steps);

	if (mode == 0)
	{
		double step = (rmax - rmin) / steps;
		
		for (uint k = 0; k <= steps; k++)
		{
			double knot = rmin + k*step;

			grid.push_back(knot);
		}
	}

	else if(mode == 1)
	{
		double emax = log(rmax+1);
		double step = emax / steps;
		
		for (uint k = 0; k <= steps; k++) 
		{
			double knot = exp(step*k)-1;
			grid.push_back(knot);
		}

	}
	else
	{
		cout << "invalid grid mode " << endl;
	}

	ofstream fs;
	fs.open(filename);
	
	for (auto& el : grid)
		fs << el << endl;
	
	fs.close();
}

double gaussleft (uint i, uint j, int l, uint o, vector <double>* absc, vector <double>* w, vector <double>* grid, bsplines* splines)
{
	double accval = 0;

	for (uint u = max(i,j); u <= min(i,j)+o-1; u++)
	{
		vector <double>* knots = splines->getknots();
		double pf = 0.5*(knots->at(u+1) - knots->at(u));

		for (uint v = 0; v < w->size(); v++)
		{
			double r = absc->at(v)*0.5*(knots->at(u+1)-knots->at(u)) + 0.5*(knots->at(u+1)+knots->at(u));
			splines->calcspline(r);
			
			if (knots->at(u) == knots->at(u+1))
			{
				accval += 0;
				break;
			}
		
			grid->push_back(r);
			
			double bi = splines->getvals(o-1,i)->at(0);
			double dj = splines->getd2(o-1,j)->at(0);
			double bj = splines->getvals(o-1,j)->at(0);

			double el = pf * w->at(v) * (bi * (-0.5*dj+(l*(l+1)/(2*pow(r,2))-1/r)*bj));

			accval += el;
		}
	}

	grid->push_back(0);

	return accval;
}


double gaussright (uint i, uint j, uint o, vector <double>* absc, vector <double>* w, bsplines* splines)
{
	double accval = 0 ;

	for (uint u = max(i,j); u <= min(i,j)+o-1; u++)
	{
		vector <double>* knots = splines->getknots() ;
		double pf = 0.5*(knots->at(u+1) - knots->at(u)) ;
		
		for (uint v = 0; v < w->size(); v++)
		{
			double r = absc->at(v)*0.5*(knots->at(u+1)-knots->at(u)) + 0.5*(knots->at(u+1)+knots->at(u)) ;
			splines->calcspline(r) ;
			
			if (knots->at(u) == knots->at(u+1))
			{
				accval += 0 ;
				break ;
			}
			
			double bi = splines->getvals(o-1,i)->at(0) ;
			double bj = splines->getvals(o-1,j)->at(0) ;
			
			double el = pf * w->at(v) *  bi * bj ;

			accval += el ;
		}
	}

	return accval ;
}

void buildmatrix (int l, uint o, vector <double>* absc, vector <double>* w, vector <double>* grid, arma::mat* lhs, arma::mat* rhs, bsplines* splines)
{
	for (uint i = 0; i < lhs->n_rows; i++)
		for (uint j = 0; j < lhs->n_cols; j++)
		{
			(*lhs)(i,j) = gaussleft(i,j,l,o,absc,w,grid,splines);
		}
	for (uint i = 0; i < rhs->n_rows; i++)
		for (uint j = 0; j < rhs->n_cols; j++)
		{
			(*rhs)(i,j) =  gaussright(i,j,o,absc,w,splines);
		}
}

void makegrid (vector <double>* absc, vector <double>* w, vector <double>* grid, vector <double>* wgrid, bsplines* splines)
{
	vector <double>* knots = splines->getknots();
	
	for (uint u = 0; u < knots->size()-1; u++)
	{
		double pf = 0.5*(knots->at(u+1) - knots->at(u));
		
		for (uint v = 0; v < w->size(); v++)
		{
			double r = pf*absc->at(v)*0.5*(knots->at(u+1)-knots->at(u)) + 0.5*(knots->at(u+1)+knots->at(u));
			splines->calcspline(r);
			
			if (knots->at(u) == knots->at(u+1))
				break;
		
			grid->push_back(r);
			wgrid->push_back(w->at(v));
		}
	}
}

double gaussleft2 (uint i, uint j, int l, vector <double>* wgrid, uint o, bsplines* splines)
{
	double accval = 0;

	vector <double>* grid = splines->getgrid(); 

	for (auto& r : (*grid))
	{
		uint k = &r - &(*grid)[0];

		double bi = splines->getvals(o-1,i)->at(k);
		double dj = splines->getd2(o-1,j)->at(k);
		double bj = splines->getvals(o-1,j)->at(k);

		double el = wgrid->at(k) * (bi * (-0.5*dj+(l*(l+1)/(2*pow(r,2))-1/r)*bj));

		accval += el;
	}

	return accval;
}

double gaussright2 (uint i, uint j, vector <double>* wgrid, uint o, bsplines* splines)
{
	double accval = 0;

	vector <double>* grid = splines->getgrid(); 

	for (auto& r : (*grid))
	{
		uint k = &r - &(*grid)[0];

		double bi = splines->getvals(o-1,i)->at(k);
		double bj = splines->getvals(o-1,j)->at(k);

		double el = wgrid->at(k) * (bi * bj);

		accval += el;
	}

	return accval;
}


void buildmatrix2 (int l, uint o, vector <double>* wgrid, arma::mat* lhs, arma::mat* rhs, bsplines* splines)
{
	for (uint i = 0; i < lhs->n_rows; i++)
		for (uint j = 0; j < lhs->n_cols; j++)
		{
			(*lhs)(i,j) = gaussleft2 (i,j,l,wgrid,o,splines);
		}
	for (uint i = 0; i < rhs->n_rows; i++)
		for (uint j = 0; j < rhs->n_cols; j++)
		{
			(*rhs)(i,j) = gaussright2 (i,j,wgrid,o,splines);
		}
}

void boundaryconds(arma::mat* lhs, arma::mat* rhs)
{
	arma::rowvec ltop(rhs->n_cols,arma::fill::zeros);
	arma::rowvec lbot(rhs->n_cols,arma::fill::zeros);
	
	(*lhs).row(0) = ltop;
	(*lhs).row(lhs->n_rows-1) = lbot;

	arma::rowvec rtop(rhs->n_cols,arma::fill::zeros);
	arma::rowvec rbot(rhs->n_cols,arma::fill::zeros);
	rtop(0) = 1;
	rbot(rbot.size()-1) = 1;

	(*rhs).row(0) = rtop;
	(*rhs).row(rhs->n_rows-1) = rbot;
}

double elvls (uint n)
{
	double E0 = -13.6 / 27.211385; // 1 Eh is 27.211385 is eV;
	double en = E0 / pow(n,2);

	return en;
}

double psi (uint n, int l, double r)
{
	double a = 1;

	double t111 = factorial(n-l-1);
	double t112 = 2*n*factorial(n+l);
	double t11 = t111/t112;
	double t1 = sqrt(pow(2/(n*a),3)*t11);
	double t2 = exp(-r/(n*a))*pow(2*r/(n*a),l);
	double t3 = pow(boost::math::laguerre(n-l-1,2*r/(n*a)),2*l+1);

	double psi = t1*t2*t3;

	//cout << "t1: " << t1 << "	t2: " << t2 << "	t3: " << t3 << endl;

	return psi;
}

void writevals(double rmin, double rmax, uint n, int l, uint i, uint steps, uint o, arma::mat* eigenvecs, arma::mat* S, bsplines* splines, string outfile)
{
	vector <double> exvals;
	exvals.reserve(steps);
	vector <double> plotgrid;
	
	double step = (rmax - rmin)/ steps;
	
	for (uint k = 0; k < steps; k++)
	{
		double r = rmin + k*step;

		double x = psi(n,l,r);

		exvals.push_back(x);
		plotgrid.push_back(r);
	}

	uint nosplines = splines->getnosplines(o-1);
	vector <double> estvals;
	estvals.reserve(steps);

	for (uint k = 0; k < steps; k++)
	{
		double r = rmin + k*step;
		
		double x = 0;

		for (uint m = 0; m < nosplines; m++)
		{	
			splines->calcspline(r);
			x += (*eigenvecs)(m,i)*splines->getvals(o-1,m)->at(0);
		}

		estvals.push_back(x);
	}

	double norm = 0;
	for(uint k = 0; k < nosplines; k++)
		for(uint j = 0; j < nosplines; j++)
		{
			norm += (*eigenvecs)(k,i)*(*eigenvecs)(j,i)*(*S)(k,j);
		}
	
	norm = sqrt(norm);	// this is the factor the norm could off by: 2^(-3/2)

	for(auto& el : estvals)
		el /= norm;

	auto sumvec = estvals;
	double csum = 0;
	for (auto& el : sumvec)
	{
		el = pow(el,2);			// square psi
		csum += el * step;		// step coult be off by factor 8? 
	}

	cout << csum << endl;

	fstream fs;
	fs.open(outfile,fstream::out);

	for (uint k = 0; k < steps; k++)
	{
		fs << exvals.at(k) << "," << plotgrid[k] << "," << pow(estvals.at(k),2) << endl;
	}

	fs.close();
}

void writegrid (vector <double>* grid, string gridfile)
{
	fstream fs;
	fs.open(gridfile,fstream::out);

	for (auto& el : (*grid))
		fs << el << endl;

	fs.close();
}

vector <double> getabsc()
{
	// could be calculated as roots to legendre polynomials, but n abscissas is sufficient for a
	// degree 2n -1 polynomial, hence we will just define five points and be fine

	double xn2 = -1.0/3.0 * sqrt(5+2*sqrt(10.0/7.0));
	double xn1 = -1.0/3.0 * sqrt(5-2*sqrt(10.0/7.0));
	double x0  =  0;
	double xp1 =  1.0/3.0 * sqrt(5-2*sqrt(10.0/7.0));
	double xp2 =  1.0/3.0 * sqrt(5+2*sqrt(10.0/7.0));

	vector <double> absc1 = {xn2,xn1,x0,xp1,xp2};
	vector <double> absc2 = {-0.973906528517172,-0.865063366688985,-0.679409568299024,-0.433395394129247,-0.148874338981631,0.148874338981631,0.433395394129247,0.679409568299024,0.865063366688985,0.973906528517172};

	return absc2;
}

vector <double> getw()
{
	// could be calculated in terms of the abscissas and derivatives of legendre polynomials, but n abscissas
	// is sufficient for a degree 2n - 1 polynomial, hence we will just define five weights and be fine

	double wn2 = (322+13*sqrt(70))/900.0;
	double wn1 = (322-13*sqrt(70))/900.0;
	double w0  = 128.0/225.0;
	double wp1 = (322-13*sqrt(70))/900.0;
	double wp2 = (322+13*sqrt(70))/900.0;

	vector <double> w1 = {wn2,wn1,w0,wp1,wp2};
	vector <double> w2 = {0.295524224714753,0.269266719309996,0.219086362515982,0.149451349150581,0.066671344308688,0.066671344308688,0.149451349150581,0.219086362515982,0.269266719309996,0.295524224714753};

	return w2;
}

int main()
{
	string knotfile = "knots.dat";
	string gridfile = "grid.dat";
	string outfile = "output.csv";

	uint order = 5;
	uint gridpts = 150;
	double tol = 1e-9;
	int l = 1;
	uint n = 3;

	double rmin = 0.0;
	double rmax = 30.0;
	uint mode = 1;

	vector <double> absc = getabsc();
	vector <double> w = getw();
	vector <double> grid;
	vector <double> wgrid;

	makeknots (rmin,rmax,gridpts,mode,knotfile);

	bsplines splines(knotfile,"",order,tol);
	uint nosplines = splines.getnosplines(order-1);
	
	makegrid (&absc,&w,&grid,&wgrid,&splines);
	writegrid(&grid,gridfile);
	bsplines splines2(knotfile,gridfile,order,tol);
	splines2.calcsplines();

	arma::mat lhs(nosplines,nosplines,arma::fill::zeros); 
	arma::mat rhs(nosplines,nosplines,arma::fill::zeros);

	buildmatrix (l,order,&absc,&w,&grid,&lhs,&rhs,&splines);
	boundaryconds (&lhs,&rhs);
	
	//cout << "lhs content:" << endl;
	//lhs.print();
	//cout << "rhs content:" << endl;
	//rhs.print();

	arma::cx_mat ceigvecs;
	arma::cx_vec ceigvals;

	eig_pair(ceigvals,ceigvecs,lhs,rhs);

	arma::vec eigvals = real(ceigvals);
	arma::mat eigvecs = real(ceigvecs);

	uint i = eigvals.size()-8;		// picking vector, should be the position of the bound state you want to study
	
	cout << "eigenvalues " << endl;
	eigvals.print();

	writevals(rmin,rmax,n,l,i,gridpts,order,&eigvecs,&rhs,&splines,outfile);
}





























