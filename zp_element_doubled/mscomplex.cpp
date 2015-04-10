#include <vector>
#include <fstream>
#include "mscomplex.h"
#include "intetration_line_tracer.h"
#include "critical_point_finder.h"
using namespace std;
namespace msc2d
{
	MSComplex2D::MSComplex2D(){}
	MSComplex2D::~MSComplex2D(){}

	bool MSComplex2D::createMSComplex2D()
	{
		CPFinder cp_finder(*this);
		cp_finder.findCriticalPoints();

		ILTracer il_tracer(*this);
		il_tracer.traceAscendingPath();
		return true;
	}
	bool MSComplex2D::saveMSComplex()const
	{
		ofstream cp("D:\\201test.txt");
		for (int i = 0; i < cp_vec.size(); ++i)
		{
			
			if (cp_vec[i].type == MAXIMAL)
			{
				cp << cp_vec[i].xy_local.first << " " << cp_vec[i].xy_local.second << " " << "MAXIMAL"<<endl;
			}
			if (cp_vec[i].type == MINIMAL)
			{
				cp << cp_vec[i].xy_local.first << " " << cp_vec[i].xy_local.second << " " << "MINIMAL" << endl;
			}
			if (cp_vec[i].type == SADDLE)
			{
				cp << cp_vec[i].xy_local.first << " " << cp_vec[i].xy_local.second << " " << "SADDLE"<<endl;
			}
			
		}
		return true;
	}
}