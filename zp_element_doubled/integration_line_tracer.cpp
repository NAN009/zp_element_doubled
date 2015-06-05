#include "integration_line_tracer.h"
#include "critical_point_finder.h"
#include <iostream>
#include <vector>
//#include "vtkReader.h"
#include "mscomplex.h"
#include<memory>
const int vr_size_1 = 1024;
const int vr_size_2 = 768;
//shared_ptr<vtkReader> vr; 
//vtkReader vr1;

using namespace std;
namespace msc2d
{
	ILTracer::ILTracer(MSComplex2D& _msc) :msc(_msc){}
	ILTracer::~ILTracer(){}

	int ILTracer::Round(double r)
	{
		return r > 0?ceil(r) :floor(r);
	//	return floor(r+0.5);
	}
	bool ILTracer::traceIntegrationPath()
	{
		cout << "Trace line path!" << endl;
		//ascend first part
		for (vector<CriticalPoint>::iterator it = msc.saddles.begin(); it != msc.saddles.end(); ++it)
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
						curr_vid = xy.first*vr_size_1 + xy.second;
					}
					if (msc.cp_vec[curr_vid].xy_local.first >= vr_size_1 - 2 || msc.cp_vec[curr_vid].xy_local.second >= vr_size_2 - 2
						|| msc.cp_vec[curr_vid].xy_local.first >= vr_size_1 - 3 || msc.cp_vec[curr_vid].xy_local.second >= vr_size_2 - 3 ||
						msc.cp_vec[curr_vid].xy_local.first <= 2 || msc.cp_vec[curr_vid].xy_local.second <= 2)
						break;
					pair<int, int> tmp_xy;
					tmp_xy = getGradDirectionUp(msc.cp_vec[curr_vid].xy_local);
					
					mesh_path.push_back(tmp_xy);
					

					prev_vid = curr_vid;
					curr_vid = tmp_xy.first*vr_size_1 + tmp_xy.second;
					if (msc.cp_vec[curr_vid].type == MINIMAL || msc.cp_vec[curr_vid + 1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1 + 1].type == MINIMAL || msc.cp_vec[curr_vid - 1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1 - 1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1 + 1].type == MINIMAL ||
						msc.cp_vec[curr_vid + vr_size_1 - 1].type == MINIMAL)
					{

						cout << "ERROR,This line is ascend！" << endl;
						break;
					}
				}
				il.startIndex = it->xy_local;
				il.endIndex = mesh_path[mesh_path.size() - 1];
			}
		}	
		//ascend second part
		for (vector<CriticalPoint>::iterator it = msc.saddles.begin(); it != msc.saddles.end(); ++it)
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

						double eig_vector_x = -msc.cp_vec[curr_vid].eig_vector1.first;
						double eig_vector_y = -msc.cp_vec[curr_vid].eig_vector1.second;
						xy = getTheSaddleBeginDirection(xy, make_pair(eig_vector_x, eig_vector_y));

						mesh_path.push_back(xy);
						curr_vid = xy.first*vr_size_1 + xy.second;
					}
					if (msc.cp_vec[curr_vid].xy_local.first >= vr_size_1 - 2 || msc.cp_vec[curr_vid].xy_local.second >= vr_size_2 - 2
						|| msc.cp_vec[curr_vid].xy_local.first >= vr_size_1 - 3 || msc.cp_vec[curr_vid].xy_local.second >= vr_size_2 - 3 ||
						msc.cp_vec[curr_vid].xy_local.first <= 2 || msc.cp_vec[curr_vid].xy_local.second <= 2)
						break;
					pair<int, int> tmp_xy;
					tmp_xy = getGradDirectionDown(msc.cp_vec[curr_vid].xy_local);
					
					mesh_path.push_back(tmp_xy);
					

					prev_vid = curr_vid;
					curr_vid = tmp_xy.first*vr_size_1 + tmp_xy.second;
					if (msc.cp_vec[curr_vid].type == MINIMAL || msc.cp_vec[curr_vid + 1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1 + 1].type == MINIMAL || msc.cp_vec[curr_vid - 1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1 - 1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1 + 1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1 - 1].type == MINIMAL)
					{

						cout << "ERROR,This line is ascend！" << endl;
						break;
					}
				}
				il.startIndex = it->xy_local;
				il.endIndex = mesh_path[mesh_path.size() - 1];
			}
		}
		//descend first part
		for (vector<CriticalPoint>::iterator it = msc.saddles.begin(); it != msc.saddles.end(); ++it)
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

						double eig_vector_x = msc.cp_vec[curr_vid].eig_vector2.first;
						double eig_vector_y = msc.cp_vec[curr_vid].eig_vector2.second;
						xy = getTheSaddleBeginDirection(xy, make_pair(eig_vector_x, eig_vector_y));

						mesh_path.push_back(xy);
						curr_vid = xy.first*vr_size_1 + xy.second;
					}
					if (msc.cp_vec[curr_vid].xy_local.first >= vr_size_1 - 2 || msc.cp_vec[curr_vid].xy_local.second >= vr_size_2 - 2
						|| msc.cp_vec[curr_vid].xy_local.first >= vr_size_1 - 3 || msc.cp_vec[curr_vid].xy_local.second >= vr_size_2 - 3 ||
						msc.cp_vec[curr_vid].xy_local.first <= 2 || msc.cp_vec[curr_vid].xy_local.second <= 2)
						break;
					pair<int, int> tmp_xy;
					tmp_xy = getGradDirectionUp1(msc.cp_vec[curr_vid].xy_local);
					
					mesh_path.push_back(tmp_xy);
					
					prev_vid = curr_vid;
					curr_vid = tmp_xy.first*vr_size_1 + tmp_xy.second;
					if (msc.cp_vec[curr_vid].type == MINIMAL || msc.cp_vec[curr_vid + 1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1 + 1].type == MINIMAL || msc.cp_vec[curr_vid - 1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1 - 1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1 + 1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1 - 1].type == MINIMAL)
					{

						cout << "This line is descend！" << endl;
						break;
					}
				}
				il.startIndex = it->xy_local;
				il.endIndex = mesh_path[mesh_path.size() - 1];
			}
		}
		//descend second part
		for (vector<CriticalPoint>::iterator it = msc.saddles.begin(); it != msc.saddles.end(); ++it)
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

						double eig_vector_x = -msc.cp_vec[curr_vid].eig_vector2.first;
						double eig_vector_y = -msc.cp_vec[curr_vid].eig_vector2.second;
						xy = getTheSaddleBeginDirection(xy, make_pair(eig_vector_x, eig_vector_y));

						mesh_path.push_back(xy);
						curr_vid = xy.first*vr_size_1 + xy.second;
					}
					if (msc.cp_vec[curr_vid].xy_local.first >= vr_size_1 - 2 || msc.cp_vec[curr_vid].xy_local.second >= vr_size_2 - 2
						|| msc.cp_vec[curr_vid].xy_local.first >= vr_size_1 - 3 || msc.cp_vec[curr_vid].xy_local.second >= vr_size_2 - 3 ||
						msc.cp_vec[curr_vid].xy_local.first <= 2 || msc.cp_vec[curr_vid].xy_local.second <= 2)
						break;
					pair<int, int> tmp_xy;
					tmp_xy = getGradDirectionDown1(msc.cp_vec[curr_vid].xy_local);
					
					mesh_path.push_back(tmp_xy);
					

					prev_vid = curr_vid;
					curr_vid = tmp_xy.first*vr_size_1 + tmp_xy.second;
					if (msc.cp_vec[curr_vid].type == MINIMAL || msc.cp_vec[curr_vid + 1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1 + 1].type == MINIMAL || msc.cp_vec[curr_vid - 1].type == MINIMAL ||
						msc.cp_vec[curr_vid - vr_size_1 - 1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1 + 1].type == MINIMAL || msc.cp_vec[curr_vid + vr_size_1 - 1].type == MINIMAL)
					{

						cout << "This line is descend！" << endl;
						break;
					}
				}
				il.startIndex = it->xy_local;
				il.endIndex = mesh_path[mesh_path.size() - 1];
			}
		}
		
		return true;
		}
	
	int Rou(double x)
	{
		if (x > 0)return 1;
		else if (x == 0)return 0;
		else return - 1;
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

		return make_pair(cur_x+eig_vector_x, cur_y+eig_vector_y);

		/*int k1x = Round(eig_vector_x);
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
		return make_pair(next_x, next_y);*/
	}
	//getGradDirection函数选择：eig_vector1选Up,-eig_vector1选down,
	//							eig_vector2选Up1 -eig_vector2选Down1
	pair<double,double> ILTracer::getGradDirectionUp(pair<int,int> xy)
	{
		//四阶龙格库塔法 h=2,待检验，Round向下取整

		int x = xy.first, y = xy.second;

		int k1x = Round(msc.cp_vec[x*vr_size_1 + y].dif.first);
		int k2x = Round(msc.cp_vec[(x + k1x)*vr_size_1 + y + 1].dif.first);
		int k3x = Round(msc.cp_vec[(x + k2x)*vr_size_1 + y + 1].dif.first);
		int k1y = Round(msc.cp_vec[x*vr_size_1 + y].dif.second);
		int k2y = Round(msc.cp_vec[(x + 1)*vr_size_1 + y + k1y].dif.second);
		int k3y = Round(msc.cp_vec[(x + 1)*vr_size_1 + y + k2y].dif.second);
		int tmp_x = Round((msc.cp_vec[x*vr_size_1 + y].dif.first +
			2 * msc.cp_vec[(x + k1x)*vr_size_1 + y + 1].dif.first +
			2 * msc.cp_vec[(x + k2x)*vr_size_1 + y + 1].dif.first +
			msc.cp_vec[(x + k3x)*vr_size_1 + y + 2].dif.first) / 3);
		int tmp_y = Round((msc.cp_vec[x*vr_size_1 + y].dif.second +
			2 * msc.cp_vec[(x + 1)*vr_size_1 + y + k1y].dif.second +
			2 * msc.cp_vec[(x + 1)*vr_size_1 + y + k2y].dif.second +
			msc.cp_vec[(x + 2)*vr_size_1 + y + k3y].dif.second) / 3);
		
		int next_X = x + tmp_x/sqrt(tmp_x*tmp_x + tmp_y*tmp_y);
		int next_Y = y + tmp_y/sqrt(tmp_x*tmp_x + tmp_y*tmp_y);

		return make_pair(next_X, next_Y);
	}
	pair<int, int> ILTracer::getGradDirectionDown(pair<int, int> xy)
	{
		//四阶龙格库塔法 h=2,待检验，Round向下取整

		int x = xy.first, y = xy.second;
		int k1x = Round(msc.cp_vec[x*vr_size_1 + y].dif.first);
		int k2x = Round(msc.cp_vec[(x - k1x)*vr_size_1 + y - 1].dif.first);
		int k3x = Round(msc.cp_vec[(x - k2x)*vr_size_1 + y - 1].dif.first);
		int k1y = Round(msc.cp_vec[x*vr_size_1 + y].dif.second);
		int k2y = Round(msc.cp_vec[(x - 1)*vr_size_1 + y - k1y].dif.second);
		int k3y = Round(msc.cp_vec[(x - 1)*vr_size_1 + y - k2y].dif.second);
		int tmp_x = Round((msc.cp_vec[x*vr_size_1 + y].dif.first +
			2 * msc.cp_vec[(x - k1x)*vr_size_1 + y - 1].dif.first +
			2 * msc.cp_vec[(x - k2x)*vr_size_1 + y - 1].dif.first +
			msc.cp_vec[(x - k3x)*vr_size_1 + y - 2].dif.first) / 3);
		int tmp_y = Round((msc.cp_vec[x*vr_size_1 + y].dif.second +
			2 * msc.cp_vec[(x - 1)*vr_size_1 + y - k1y].dif.second +
			2 * msc.cp_vec[(x - 1)*vr_size_1 + y - k2y].dif.second +
			msc.cp_vec[(x - 2)*vr_size_1 + y - k3y].dif.second) / 3);
		
		int next_X = x + tmp_x;
		int next_Y = y + tmp_y;

		return make_pair(next_X, next_Y);

	}
	
	pair<double, double> ILTracer::getGradDirectionUp1(pair<int, int> xy)
	{
		//四阶龙格库塔法 h=2,待检验，Round向下取整

		int x = xy.first, y = xy.second;
		int k1x = Round(msc.cp_vec[x*vr_size_1 + y].dif.first);
		int k2x = Round(msc.cp_vec[(x + k1x)*vr_size_1 + y + 1].dif.first);
		int k3x = Round(msc.cp_vec[(x + k2x)*vr_size_1 + y + 1].dif.first);
		int k1y = Round(msc.cp_vec[x*vr_size_1 + y].dif.second);
		int k2y = Round(msc.cp_vec[(x + 1)*vr_size_1 + y + k1y].dif.second);
		int k3y = Round(msc.cp_vec[(x + 1)*vr_size_1 + y + k2y].dif.second);
		int tmp_x = Round((msc.cp_vec[x*vr_size_1 + y].dif.first +
			2 * msc.cp_vec[(x + k1x)*vr_size_1 + y + 1].dif.first +
			2 * msc.cp_vec[(x + k2x)*vr_size_1 + y + 1].dif.first +
			msc.cp_vec[(x + k3x)*vr_size_1 + y + 2].dif.first) / 3);
		int tmp_y = Round((msc.cp_vec[x*vr_size_1 + y].dif.second +
			2 * msc.cp_vec[(x + 1)*vr_size_1 + y + k1y].dif.second +
			2 * msc.cp_vec[(x + 1)*vr_size_1 + y + k2y].dif.second +
			msc.cp_vec[(x + 2)*vr_size_1 + y + k3y].dif.second) / 3);
		
		int next_X = x - tmp_x;
		int next_Y = y - tmp_y;

		return make_pair(next_X, next_Y);
		
	}
	pair<int, int> ILTracer::getGradDirectionDown1(pair<int, int> xy)
	{
		//四阶龙格库塔法 h=2,待检验，Round向下取整

		int x = xy.first, y = xy.second;
		int k1x = Round(msc.cp_vec[x*vr_size_1 + y].dif.first);
		int k2x = Round(msc.cp_vec[(x - k1x)*vr_size_1 + y - 1].dif.first);
		int k3x = Round(msc.cp_vec[(x - k2x)*vr_size_1 + y - 1].dif.first);
		int k1y = Round(msc.cp_vec[x*vr_size_1 + y].dif.second);
		int k2y = Round(msc.cp_vec[(x - 1)*vr_size_1 + y - k1y].dif.second);
		int k3y = Round(msc.cp_vec[(x - 1)*vr_size_1 + y - k2y].dif.second);
		int tmp_x = Round((msc.cp_vec[x*vr_size_1 + y].dif.first +
			2 * msc.cp_vec[(x - k1x)*vr_size_1 + y - 1].dif.first +
			2 * msc.cp_vec[(x - k2x)*vr_size_1 + y - 1].dif.first +
			msc.cp_vec[(x - k3x)*vr_size_1 + y - 2].dif.first) / 3);
		int tmp_y = Round((msc.cp_vec[x*vr_size_1 + y].dif.second +
			2 * msc.cp_vec[(x - 1)*vr_size_1 + y - k1y].dif.second +
			2 * msc.cp_vec[(x - 1)*vr_size_1 + y - k2y].dif.second +
			msc.cp_vec[(x - 2)*vr_size_1 + y - k3y].dif.second) / 3);

		int next_X = x - tmp_x;
		int next_Y = y - tmp_y;

		return make_pair(next_X, next_Y);

	}

}

