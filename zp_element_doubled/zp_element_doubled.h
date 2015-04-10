#pragma once

#include <vector>
#include <string>
#include <string.h>

//#include "point_evaluator.h"
#include "global_functions.h"

typedef std::vector<std::vector<std::vector<double> > > image_3d;

//class reconstructor
//{
//public:
//	reconstructor(char*,char*);
//	~reconstructor(void);
//	void set_dimension(int,int);
//	void set_step(double,double);
//	void initialize();
//	void read_data();
//	void reconstruct();
//
//	bool with_header;
//private:
//	char* original_data_path;
//	char* output_data_path;
//	//point_evaluator *pe;
//	global_functions *gf;
//	int L;
//	int W;
//	double stepx;
//	double stepy;
//	int nx;
//	int ny;
//
//	std::vector<double> x_ordinates;
//	std::vector<double> y_ordinates;
//
//	std::vector<std::vector<std::vector<double> > > data;
////	std::vector<std::vector<double> > reconstructed_piece;
//
//	
//};
