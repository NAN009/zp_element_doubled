#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include "global_functions.h"

global_functions::global_functions(void){}

global_functions::~global_functions(void){}

//double to string function.
string global_functions::double2str(double d)
{
    stringstream ss;
    ss << setprecision(16) <<d;
    return ss.str();
}

//float to string function.
string global_functions::double2str(float d)
{
    stringstream ss;
    ss << setprecision(16) <<d;
    return ss.str();
}

//integer to string function.
string global_functions::int2str(int a)
{
    stringstream ss;
    ss<<a;
    return ss.str();
}

//string to double function.
double global_functions::str2double(string s)//数据转换，输入s  输出d
{
    stringstream ss;
    double d;
    ss<<s;
    ss>>d;

    return d;
}

//string to float function.
float global_functions::str2float(string s)
{
    stringstream ss;
    float d;
    ss<<s;
    ss>>d;

    return d;
}

//string to integer function.
int global_functions::str2int(string s)
{
    stringstream ss;
    int d;
    ss<<s;
    ss>>d;

     return d;
}

//string to long integer function.
long global_functions::str2long(string s)
{
    stringstream ss;
    long d;
    ss<<s;
    ss>>d;

    return d;
}

//scientific format to not scientific format
string global_functions::sci2nsci(double s)
{
   stringstream ss;
   ss<<setprecision(16)<<fixed<<s;
   return ss.str();
}
