#pragma once

#include "point_evaluator.h"
#include "point_evaluator_7direction.h"
#include "point_evaluator_6direction1.h"
#include "point_evaluator_6direction2.h"

class box_spline
{
public:
	box_spline(void);
	~box_spline(void);
	double compute_value(double x[]);
	double compute_gradient_x(double x[]);
	double compute_gradient_y(double x[]);
	double compute_gradient_xx(double x[]);
	double compute_gradient_xy(double x[]);
	double compute_gradient_yy(double x[]);
private:
	point_evaluator* pe;
	point_evaluator_7direction* pe7;
	point_evaluator_6direction1* pe61;
	point_evaluator_6direction2* pe62;
};
