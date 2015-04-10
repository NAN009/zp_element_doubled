#pragma once

#include "global_functions.h"

class point_evaluator_7direction
{
public:
	point_evaluator_7direction(void);
	~point_evaluator_7direction(void);
	void initialize(void);
	double evaluate(double*);
private:
	global_functions *gf;
	double inv_V[3][3][4];
	double cp[2][21][4][5][6];
    double Ne[2][2];
	double bary[3];
	double x[2];

	int factorial(int);
	bool out_boundary(double*);
	double my_pow(double,int);
};
