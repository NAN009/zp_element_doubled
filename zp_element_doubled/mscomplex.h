#ifndef MSCOMPLEX_H
#define MSCOMPLEX_H
#include <vector>
#include <utility>

namespace msc2d
{using namespace  std;
	enum CriticalPointType
	{
		MINIMAL = 0,
		SADDLE,

		MAXIMAL,
		REGULAR
	};
	struct CriticalPoint
	{
		CriticalPointType type;		//0 --minimal, 1 --saddle, 2 --maximal
		pair<int, int> xy_local;				//index in original mesh
		int meshIndex;
		//缺少源文件中CriticalPointNeighborArray类型
		pair<double, double> eig_vector1;//特征向量
		pair<double, double> eig_vector2;

		pair<double, double> dif;//一阶导数<dif_x,dif_y>
	};
	typedef std::vector<CriticalPoint> CriticalPointArray;

	typedef std::vector< pair<int,int> > PATH;
	struct IntegrationLine
	{
		pair<int,int> startIndex, endIndex;//起始点x,y的位置
		PATH path;//the index into original mesh, path always start from a saddle point to a max/min point
	};
	typedef std::vector<IntegrationLine> IntegrationLineArray;
	class MSComplex2D
	{
	public:
		MSComplex2D();
		~MSComplex2D();
		
		bool createMSComplex2D();
		bool saveMSComplex()const;
	private:


		CriticalPointArray cp_vec;
		CriticalPointArray minPoints;
		CriticalPointArray maxPoints;
		CriticalPointArray saddles;
		CriticalPointArray keyPoint;

		IntegrationLineArray il_vec;

		friend class CPFinder;
		friend class ILTracer;
	};
}
#endif
