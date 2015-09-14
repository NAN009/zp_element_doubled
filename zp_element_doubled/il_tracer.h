#ifndef il_tracer_
#define il_tracer_
#include <utility>
#include <vector>
#include <iostream>
#include "mscomplex.h"
namespace msc2d
{
	class CriticalPoint;
	class ILTracer
	{
		class WEdge
		{
		public:
			//存放有所有最大值的序数
			std::vector<std::pair<size_t, size_t>> max_ranges;
			std::vector<std::pair<size_t, size_t>> min_ranges;
		};
	public:
		ILTracer(MSComplex2D& _msc);
		~ILTracer();

		bool traceIntegrationPath_RungeKutta2d();
		bool traceIntegrationPath_RungeKutta5();
		bool traceIntegrationPath();
		bool isMaximal(double x, double y);
		bool isMinimal(double x, double y);
		bool isSaddle(double x, double y);
		bool isKeyPoint(double x, double y);
		pair<double, double> ILTracer::SaddlePosition(double x, double y);
		bool isBoundary(double x, double y);
		pair<double, double> getTheSaddleBeginDirection(pair<double, double> xy, pair<double, double> eig_vector);
		pair<double, double> getGradDirectionUp(pair<double,double> xy);
		pair<double, double> getGradDirectionDown(pair<double, double> xy);
		
		int Round(double r);
	private:
		
		MSComplex2D& msc;
	};
}
#endif