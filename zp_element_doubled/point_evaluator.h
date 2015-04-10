#pragma once

#include "global_functions.h"

class point_evaluator
{
public:
	point_evaluator(void);
	~point_evaluator(void);
	void initialize(void);
	double evaluate(double*);
private:
	global_functions *gf;
	double inv_V[3][3][4];
	double cp[2][28][4][6][6];
	double p[28][4][6][6];
    double Ne[2][2];
	double bary[3];
	double x[2];

	int factorial(int);
	bool out_boundary(double*);
	double my_pow(double,int);
};
