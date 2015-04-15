#include <iostream>
#include "box_spline.h"

box_spline::box_spline(void)
{
	pe = new point_evaluator();
	pe7 = new point_evaluator_7direction();
	pe61 = new point_evaluator_6direction1();
	pe62 = new point_evaluator_6direction2();
	pe->initialize();
	pe7->initialize();
	pe61->initialize();
	pe62->initialize();
}

box_spline::~box_spline(void)
{
	delete pe;
	delete pe7;
	delete pe61;
	delete pe62;
}

double box_spline::compute_value(double x[])
{
	double* p = new double[2];
	double* tmp = p;

	*tmp = x[0];
	tmp++;
	*tmp = x[1];

	return pe->evaluate(p);
}

double box_spline::compute_gradient_x(double x[])
{
	//double* p1 = new double[2];
	//double* p2 = new double[2];
	//double* tmp1 = p1;
	//double* tmp2 = p2;
	//
	//*tmp1 = x[0]+0.5;
	//*tmp2 = x[0]-0.5;
	//tmp1++;
	//tmp2++;
	//*tmp1 = x[1];
	//*tmp2 = x[1];
	double p1[2],p2[2];
	p1[0] = x[0] + 0.5; p2[0] = x[0] - 0.5;
	p1[1] = x[1]; p2[1] = x[1];
	return pe7->evaluate(p1)-pe7->evaluate(p2);//401sin£ºpe7->evaluate(p2)-pe7->evaluate(p1)
	
	}

double box_spline::compute_gradient_y(double x[])
{
	double y[2];

	y[0] = x[1];
	y[1] = x[0];
	return this->compute_gradient_x(y);
}

double box_spline::compute_gradient_xx(double x[])
{
	double x1[2];
	double x2[2];
	double x3[2];

	x1[0] = x[0]+1;
	x1[1] = x[1];

	x2[0] = x[0]-1;
	x2[1] = x[1];

	x3[0] = x[0];
	x3[1] = x[1];

	return (pe61->evaluate(x1))+(pe61->evaluate(x2))-2*(pe61->evaluate(x3));
}

double box_spline::compute_gradient_xy(double x[])
{
	double x1[2];
	double x2[2];
	double x3[2];
	double x4[2];

	x1[0] = x[0]+0.5;
	x1[1] = x[1]+0.5;

	x2[0] = x[0]-0.5;
	x2[1] = x[1]+0.5;

	x3[0] = x[0]+0.5;
	x3[1] = x[1]-0.5;

	x4[0] = x[0]-0.5;
	x4[1] = x[1]-0.5;

	return (pe62->evaluate(x1))-(pe62->evaluate(x2))-(pe62->evaluate(x3))+(pe62->evaluate(x4));
}

double box_spline::compute_gradient_yy(double x[])
{
	double y[2];

	y[0] = x[1];
	y[1] = x[0];
	return this->compute_gradient_xx(y);
}