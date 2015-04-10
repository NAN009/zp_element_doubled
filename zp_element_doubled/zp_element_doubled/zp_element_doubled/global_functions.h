#pragma once

#include <string.h>
#include <string>

using namespace std;

class global_functions
{
public:
	global_functions(void);
	~global_functions(void);
	string double2str(double);
	string double2str(float);
	string int2str(int);
	double str2double(string);
    float str2float(string);
	int str2int(string);
	long str2long(string);
	string sci2nsci(double);
};
