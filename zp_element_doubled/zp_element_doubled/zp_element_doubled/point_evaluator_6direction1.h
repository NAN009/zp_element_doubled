#pragma once

#include "global_functions.h"

class point_evaluator_6direction1
{
public:
	point_evaluator_6direction1(void);
	~point_evaluator_6direction1(void);
    void initialize(void);
	double evaluate(double*);
private:
	global_functions *gf;
	double inv_V[3][3][4];
	double cp[2][15][4][4][6];
    double Ne[2][2];
	double bary[3];
	double x[2];

	int factorial(int);
	bool out_boundary(double*);
	double my_pow(double,int);
};
