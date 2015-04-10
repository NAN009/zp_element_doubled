#include "intetration_line_tracer.h"
#include "critical_point_finder.h"
#include <iostream>
#include <vector>
//#include "vtkReader.h"
#include "mscomplex.h"
#include<memory>
const int vr_size = 401;
//shared_ptr<vtkReader> vr; 
//vtkReader vr1;


using namespace std;
namespace msc2d
{
	ILTracer::ILTracer(MSComplex2D& _msc) :msc(_msc){}
	ILTracer::~ILTracer(){}

	
	int ILTracer::Round(double r)
	{
		if (r > 0)
			return ceil(r);
			else
			return floor(r);
		//return ceil(r);
	}
	bool ILTracer::traceAscendingPath()
	{
		cout << "Trace ascending path" << endl;
		for (vector<CriticalPoint>::iterator it = msc.cp_vec.begin(); it != msc.cp_vec.end(); ++it)
		{

			if (it->type == SADDLE)
			{
				CriticalPoint& sad = *it;

				int prev_vid = -1, curr_vid = sad.meshIndex;

				msc.il_vec.push_back(IntegrationLine());
				IntegrationLine &il = msc.il_vec[msc.il_vec.size() - 1];
				PATH& mesh_path = il.path;
				mesh_path.push_back(make_pair(sad.xy_local.first, sad.xy_local.second));

				while (msc.cp_vec[curr_vid].type != MAXIMAL)
				{
					while (msc.cp_vec[curr_vid].type == SADDLE)
					{
						pair<int, int> xy = msc.cp_vec[curr_vid].xy_local;

						double eig_vector_x = msc.cp_vec[curr_vid].eig_vector1.first;
						double eig_vector_y = msc.cp_vec[curr_vid].eig_vector1.second;
						xy = getTheSaddleBeginDirection(xy, make_pair(eig_vector_x, eig_vector_y));
					
						mesh_path.push_back(xy);
						curr_vid = xy.first*vr_size + xy.second;
					}
					pair<int, int> tmp_xy;
					tmp_xy = getGradDirection(msc.cp_vec[curr_vid].xy_local);
					if (tmp_xy.first == msc.cp_vec[curr_vid].xy_local.first&&
						tmp_xy.second == msc.cp_vec[curr_vid].xy_local.second)
					{

						cout << "ERROR,This line is ascend！" << endl;
						break;
					}
					mesh_path.push_back(tmp_xy);
					if (tmp_xy.first >= vr_size - 2 || tmp_xy.second >= vr_size - 2
						|| tmp_xy.first >= vr_size - 3 || tmp_xy.second >= vr_size - 3 ||
						tmp_xy.first <= 0 || tmp_xy.second <= 0)
						break;

					prev_vid = curr_vid;
					curr_vid = tmp_xy.first*vr_size + tmp_xy.second;
					if (msc.cp_vec[curr_vid].type == MINIMAL || msc.cp_vec[curr_vid + 1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size].type == MINIMAL || msc.cp_vec[curr_vid + vr_size].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size + 1].type == MINIMAL || msc.cp_vec[curr_vid - 1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size - 1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size + 1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size - 1].type == MINIMAL)
					{

						cout << "ERROR,This line is ascend！" << endl;
						break;
					}
				}
				il.startIndex = it->xy_local;
				il.endIndex = mesh_path[mesh_path.size() - 1];
			}
		}		
		return true;
		}
	
	pair<double,double>ILTracer::getTheSaddleBeginDirection(pair<int,int> xy,pair<double,double> eig_vector)
	{
		/*int x = xy.first, y = xy.second;

		double k1x = eig_vector.first;
		double k2x = Round(EvaluateGradX(x + k1x, y + 1));
		double k3x = Round(EvaluateGradX(x + k2x, y + 1));
		double k4x = Round(EvaluateGradX(x + k3x, y + 2));

		double k1y = eig_vector.second;
		double k2y = Round(EvaluateGradY(x + 1, y + k1y));
		double k3y = Round(EvaluateGradY(x + 1, y + k2y));
		double k4y = Round(EvaluateGradY(x + 2, y + k3y));

		double kx = x + (k1x + 2 * k2x + 2 * k3x + k4x) / 3;
		double ky = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 3;

		return make_pair(Round(kx),Round(ky));*/

		int cur_x = xy.first;
		int cur_y = xy.second;
		double eig_vector_x = eig_vector.first;
		double eig_vector_y = eig_vector.second;
		int k1x = Round(eig_vector_x);
		int k2x = Round(msc.cp_vec[(cur_x + k1x)*vr_size + cur_y + 1].dif.first);
		int k3x = Round(msc.cp_vec[(cur_x + k2x)*vr_size + cur_y + 1].dif.first);

		int k1y = Round(eig_vector_y);
		int k2y = Round(msc.cp_vec[cur_x + 1 + (cur_y + k1y)*vr_size].dif.second);
		int k3y = Round(msc.cp_vec[cur_x + 1 + (cur_y + k2y)*vr_size].dif.second);

		int next_x = cur_x + Round(
		(eig_vector_x +
		2 * (msc.cp_vec[(cur_x + k1x)*vr_size + cur_y + 1].dif.first) +
		2 * (msc.cp_vec[(cur_x + k2x)*vr_size + cur_y + 1].dif.first) +
		msc.cp_vec[(cur_x + k3x)*vr_size + cur_y + 2].dif.first) / 3);
		int next_y = cur_y + Round(
		(eig_vector_y +
		2 * (msc.cp_vec[cur_x + 1 + (cur_y + k1y)*vr_size].dif.second) +
		2 * (msc.cp_vec[cur_x + 1 + (cur_y + k2y)*vr_size].dif.second) +
		msc.cp_vec[cur_x + 2 + (cur_y + k3y)*vr_size].dif.second) / 3);
		return make_pair(next_x, next_y);
	}
	pair<double,double> ILTracer::getGradDirection(pair<int,int> xy)
	{
		//四阶龙格库塔法 h=2,待检验，Round向下取整

		int x = xy.first, y = xy.second;
/*
		double k1x = Round(EvaluateGradX(x, y));
		double k2x = Round(EvaluateGradX(x + k1x, y + 1));
		double k3x = Round(EvaluateGradX(x + k2x, y + 1));
		double k4x = Round(EvaluateGradX(x + k3x, y + 2));

		double k1y = Round(EvaluateGradY(x, y));
		double k2y = Round(EvaluateGradY(x + 1, y + k1y));
		double k3y = Round(EvaluateGradY(x + 1, y + k2y));
		double k4y = Round(EvaluateGradY(x + 2, y + k3y));

		double kx = x + (k1x + 2*k2x + 2*k3x + k4x) / 3;
		double ky = y + (k1y + 2*k2y + 2*k3y + k4y) / 3;

		return make_pair(Round(kx), Round(ky));*/

		int k1x = Round(msc.cp_vec[x*vr_size + y].dif.first);
		int k2x = Round(msc.cp_vec[(x + k1x)*vr_size + y + 1].dif.first);
		int k3x = Round(msc.cp_vec[(x + k2x)*vr_size + y + 1].dif.first);
		int k1y = Round(msc.cp_vec[y*vr_size + x].dif.second);
		int k2y = Round(msc.cp_vec[(y + k1y)*vr_size + x + 1].dif.second);
		int k3y = Round(msc.cp_vec[(y + k2y)*vr_size + x + 1].dif.second);
		int tmp_x = Round((msc.cp_vec[x*vr_size + y].dif.first +
			2 * msc.cp_vec[(x + k1x)*vr_size + y + 1].dif.first +
			2 * msc.cp_vec[(x + k2x)*vr_size + y + 1].dif.first +
			msc.cp_vec[(x + k3x)*vr_size + y + 2].dif.first) / 3);
		int tmp_y = Round((msc.cp_vec[y*vr_size + x].dif.second +
			2 * msc.cp_vec[(y + k1y)*vr_size + x + 1].dif.second +
			2 * msc.cp_vec[(y + k2y)*vr_size + x + 1].dif.second +
			msc.cp_vec[(y + k3y)*vr_size + x + 2].dif.second) / 3);
		if (tmp_x < 0 || tmp_y < 0)
			return make_pair(x, y);
		int next_X = x + tmp_x;
		int next_Y = y + tmp_y;

		return make_pair(next_X, next_Y);

	}
	pair<int, int> ILTracer::getGradDirection1(pair<int, int> xy)
	{
		//四阶龙格库塔法 h=2,待检验，Round向下取整

		int x = xy.first, y = xy.second;
		int k1x = Round(msc.cp_vec[x*vr_size + y].dif.first);
		int k2x = Round(msc.cp_vec[(x - k1x)*vr_size + y - 1].dif.first);
		int k3x = Round(msc.cp_vec[(x - k2x)*vr_size + y - 1].dif.first);
		int next_X = x + Round((msc.cp_vec[x*vr_size + y].dif.first +
		2 * msc.cp_vec[(x - k1x)*vr_size + y - 1].dif.first +
		2 * msc.cp_vec[(x - k2x)*vr_size + y - 1].dif.first +
		msc.cp_vec[(x - k3x)*vr_size + y - 2].dif.first) / 3);
		int k1y = Round(msc.cp_vec[y*vr_size + x].dif.second);
		int k2y = Round(msc.cp_vec[(y - k1y)*vr_size + x - 1].dif.second);
		int k3y = Round(msc.cp_vec[(y - k2y)*vr_size + x - 1].dif.second);
		int next_Y = y + Round((msc.cp_vec[y*vr_size + x].dif.second +
		2 * msc.cp_vec[(y - k1y)*vr_size + x - 1].dif.second +
		2 * msc.cp_vec[(y - k2y)*vr_size + x - 1].dif.second +
		msc.cp_vec[(y - k3y)*vr_size + x - 2].dif.second) / 3);

		return make_pair(next_X, next_Y);

	}
}

