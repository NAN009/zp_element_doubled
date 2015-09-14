#include <vector>
#include <fstream>
#include "mscomplex.h"
#include <vector>
#include "il_tracer.h"
//#include "integration_line_tracer.h"
//#include "critical_point_finder.h"
#include "criticalPoint.h"
using namespace std;
namespace msc2d
{
	MSComplex2D::MSComplex2D(){}
	MSComplex2D::~MSComplex2D(){}

	bool MSComplex2D::createMSComplex2D()
	{
		//CPFinder cp_finder(*this);
		//cp_finder.findCriticalPoints();
		CriticalPointFind cp_finder(*this);
		cp_finder.findCriticalPoint();

		ILTracer il_tracer(*this);
		il_tracer.traceIntegrationPath_RungeKutta5();
		//il_tracer.traceIntegrationPath();
		return true;
	}
	bool MSComplex2D::saveMSComplex()const
	{
		ofstream il("D:\\newData\\taiwan_2\\taiwan_2_gray_filter_il.txt");
		for (int i = 0; i < il_vec.size(); ++i)
		{
			//il<<il_vec[i].startIndex.first << " " << il_vec[i].endIndex;
			for (int j = 0; j < il_vec[i].path.size();++j)
			{
				il << il_vec[i].path[j].first << " " << il_vec[i].path[j].second << endl;
			}
			il << endl;
		}
		return true;
	}
}