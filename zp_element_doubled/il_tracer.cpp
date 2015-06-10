#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "criticalPoint.h"
#include "il_tracer.h"
#include "mscomplex.h"
#include<memory>
#define  MAXSTEP  1
#define EPS 1e-3
#define NORM2(x,y) ((x)*(x)+(y)*(y))
inline double norm(double x, double y)
{
	return x*x + y*y;
}
const int vr_size_1 = 214;
const int vr_size_2 = 286;
const double LARGE_ZERO_EPSILON = 0.1;
const double APPROACH = 0.1;
using namespace std;
namespace msc2d
{
	ILTracer::ILTracer(MSComplex2D& _msc):msc(_msc){}
	ILTracer::~ILTracer(){}
	bool ILTracer::traceIntegrationPath()
	{
		cout << "Trace line path!" << endl;
		for (vector<CriticalPoint>::iterator it = msc.saddles.begin(); it != msc.saddles.end(); ++it)
		{
			double px, py, posx, posy, dstepx, dstepy, testx, testy;
			double pval, testval, step, l, dist, dir_x, dir_y;
			unsigned int path, i, j, num;
			for (path = 0; path < 2; ++path)
			{
				msc.il_vec.push_back(IntegrationLine());
				IntegrationLine &il = msc.il_vec[msc.il_vec.size() - 1];
				PATH& mesh_path = il.path;
				num = 0;
				posx = it->xy_local.first; posy = it->xy_local.second;
				mesh_path.push_back(make_pair(posx, posy));
				px = it->xy_local.first; py = it->xy_local.second;
				testx = it->xy_local.first; testy = it->xy_local.second;

				dir_x = it->eig_vector1.first;
				dir_y = it->eig_vector1.second;
				if (path > 0)
				{
					dir_x = -dir_x;
					dir_y = -dir_y;
				}
				while (1)
				{
					posx = px + dir_x*MAXSTEP;
					posy = py + dir_y*MAXSTEP;
					if (isBoundary(posx,posy))						
					{
						il.startIndex = it->xy_local;
						il.endIndex = mesh_path[mesh_path.size() - 1];
						break;
					}
					pval = EvaluateVal(px, py);
					step = MAXSTEP;
					while (step >= EPS)
					{
						posx = px + dir_x*step;
						posy = py + dir_y*step;
						testx = posx; testy = posy;
						testval = EvaluateVal(testx, testy);
						step *= 0.618;
						if (testval > pval)
						{
							posx = px + dir_x*step;
							posy = py + dir_y*step;
							px = posx; py = posy;
							if (pval >= EvaluateVal(px, py))
								px = testx; py = testy;
							dir_x = EvaluateGradX(px, py);
							dir_y = EvaluateGradY(px, py);
							//l = NORM2(dir_x, dir_y);
							//l = norm(dir_x, dir_y);
							//if (l > 0)
							//{
							//	//l = 1.0 /sqrt(l);									
							//	dir_x *= l;
							//	dir_y *= l;
							//}
							break;
						}
					}
					if (step < EPS)
					{
						dist = MAXSTEP;
						for (int i = 0; i < msc.maxPoints.size(); ++i)
						{
							dstepx = msc.maxPoints[i].xy_local.first - px;
							dstepy = msc.maxPoints[i].xy_local.second - py;
							if (norm(dstepx, dstepy) < dist)
							{
								j = i;
								dist = norm(dstepx, dstepy);
							}
						}
						if (dist < MAXSTEP)
						{
							posx = msc.maxPoints[j].xy_local.first;
							posy = msc.maxPoints[j].xy_local.second;
							mesh_path.push_back(make_pair(posx, posy));
							il.startIndex = it->xy_local;
							il.endIndex = mesh_path[mesh_path.size() - 1];
						}
						else
						{
							CriticalPoint newmax;
							newmax.xy_local = make_pair(posx, posy);
							newmax.type = CriticalPointType::MAXIMAL;
							msc.maxPoints.push_back(newmax);
							msc.keyPoint.push_back(newmax);
							mesh_path.push_back(make_pair(posx, posy));
							il.startIndex = it->xy_local;
							il.endIndex = mesh_path[mesh_path.size() - 1];
						}
						break;
					}
					else
					{
						dstepx = mesh_path[num].first - posx;
						dstepy = mesh_path[num].second - py;
						if (norm(dstepx, dstepy) > MAXSTEP*MAXSTEP)
						{
							mesh_path.push_back(make_pair(posx, posy));
							++num;
						}
					}
				}
			}
		}
		return true;
	}
	bool ILTracer::traceIntegrationPath_RungeKutta()
	{
		ofstream il_out("D:\\newData\\2Ddata\\gd1wie120_20_iltext.txt");
		cout << "Trace line path!" << endl;	
	
		for (vector<CriticalPoint>::iterator it = msc.saddles.begin(); it != msc.saddles.end(); ++it)
		{
			if (it->type == SADDLE)
			{
				CriticalPoint& sad = *it;
				pair<double, double> curr(sad.xy_local);
				int prev_vid = -1, curr_vid = sad.meshIndex;
				msc.il_vec.push_back(IntegrationLine());
				IntegrationLine &il = msc.il_vec[msc.il_vec.size() - 1];
				PATH& mesh_path = il.path;
				mesh_path.push_back(sad.xy_local);
				il_out << sad.xy_local.first << " " << sad.xy_local.second << endl;
				while (!isMaximal(curr.first, curr.second))
				{
					if (isBoundary(curr.first, curr.second))
						break;
					while (isSaddle(curr.first, curr.second))
					{
						pair<double, double> saddle_positon = SaddlePosition(curr.first, curr.second);
						pair<double, double> xy = curr;
						double eig_vector_x = msc.cp_vec[saddle_positon.first*vr_size_2 + saddle_positon.second].eig_vector1.first;
						double eig_vector_y = msc.cp_vec[saddle_positon.first*vr_size_2 + saddle_positon.second].eig_vector1.second;
						xy = getTheSaddleBeginDirection(xy, make_pair(eig_vector_x, eig_vector_y));
						mesh_path.push_back(xy);
						il_out << xy.first << " " << xy.second << endl;
						//curr_vid = Round(xy.first)*vr_size_2 + Round(xy.second);
						curr = xy;
					}
					if (isBoundary(curr.first, curr.second))
						break;
					pair<double, double> tmp_xy;
					tmp_xy = getGradDirectionDown(curr);

					if (norm(tmp_xy.first - mesh_path[mesh_path.size() - 2].first, tmp_xy.second - mesh_path[mesh_path.size() - 2].second) < APPROACH*APPROACH || norm(tmp_xy.first - mesh_path[mesh_path.size() - 1].first, tmp_xy.second - mesh_path[mesh_path.size() - 1].second) < APPROACH*APPROACH)
					{
						if ((EvaluateGradX(tmp_xy.first, tmp_xy.second) < LARGE_ZERO_EPSILON&&
							EvaluateGradY(tmp_xy.first, tmp_xy.second) < LARGE_ZERO_EPSILON) ||
							isKeyPoint(tmp_xy.first, tmp_xy.second))
						{
							mesh_path.push_back(tmp_xy);
							il_out << tmp_xy.first << " " << tmp_xy.second << endl;
							//prev_vid = curr_vid;
							//curr_vid = Round(tmp_xy.first)*vr_size_2 + Round(tmp_xy.second);
							curr = tmp_xy;
							cout << "This endPoint is a CriticalPoint！" << endl;
							break;
						}
						else
						{
							double x_dir = EvaluateGradX(tmp_xy.first, tmp_xy.second);
							double y_dir = EvaluateGradY(tmp_xy.first, tmp_xy.second);
							tmp_xy.first += x_dir / sqrt(norm(x_dir, y_dir));
							tmp_xy.second += y_dir / sqrt(norm(x_dir, y_dir));
							bool b = true;
							for (int t = mesh_path.size() - 4; t < mesh_path.size(); ++t)
							{
								if (tmp_xy.first - mesh_path[t].first>1 || tmp_xy.second - mesh_path[t].second > 1)
								{
									b = false;
								}
							}
							if (b == true)
								break;
						}
					}
					mesh_path.push_back(tmp_xy);
					il_out << tmp_xy.first << " " << tmp_xy.second << endl;
					//prev_vid = curr_vid;
					//curr_vid = Round(tmp_xy.first)*vr_size_2 + Round(tmp_xy.second);
					curr = tmp_xy;
					if (isMinimal(curr.first, curr.second))
					{
						cout << "ERROR,This line is ascend！" << endl;
						break;
					}
				}
				il.startIndex = it->xy_local;
				il.endIndex = mesh_path[mesh_path.size() - 1];
				il_out << endl;
			}
		}
		for (vector<CriticalPoint>::iterator it = msc.saddles.begin(); it != msc.saddles.end(); ++it)
		{
			if (it->type == SADDLE)
			{
				CriticalPoint& sad = *it;
				pair<double, double> curr(sad.xy_local);
				int prev_vid = -1, curr_vid = sad.meshIndex;
				msc.il_vec.push_back(IntegrationLine());
				IntegrationLine &il = msc.il_vec[msc.il_vec.size() - 1];
				PATH& mesh_path = il.path;
				mesh_path.push_back(sad.xy_local);
				il_out << sad.xy_local.first << " " << sad.xy_local.second << endl;
				while (!isMaximal(curr.first, curr.second))
				{
					if (isBoundary(curr.first, curr.second))
						break;
					while (isSaddle(curr.first, curr.second))
					{
						pair<double, double> saddle_positon = SaddlePosition(curr.first, curr.second);
						pair<double, double> xy = curr;
						double eig_vector_x = msc.cp_vec[saddle_positon.first*vr_size_2 + saddle_positon.second].eig_vector2.first;
						double eig_vector_y = msc.cp_vec[saddle_positon.first*vr_size_2 + saddle_positon.second].eig_vector2.second;
						xy = getTheSaddleBeginDirection(xy, make_pair(eig_vector_x, eig_vector_y));
						mesh_path.push_back(xy);
						il_out << xy.first << " " << xy.second << endl;
						//curr_vid = Round(xy.first)*vr_size_2 + Round(xy.second);
						curr = xy;
					}
					if (isBoundary(curr.first, curr.second))
						break;
					pair<double, double> tmp_xy;
					tmp_xy = getGradDirectionDown(curr);

					if (norm(tmp_xy.first - mesh_path[mesh_path.size() - 2].first, tmp_xy.second - mesh_path[mesh_path.size() - 2].second) < APPROACH*APPROACH || norm(tmp_xy.first - mesh_path[mesh_path.size() - 1].first, tmp_xy.second - mesh_path[mesh_path.size() - 1].second) < APPROACH*APPROACH)
					{
						if ((EvaluateGradX(tmp_xy.first, tmp_xy.second) < LARGE_ZERO_EPSILON&&
							EvaluateGradY(tmp_xy.first, tmp_xy.second) < LARGE_ZERO_EPSILON) ||
							isKeyPoint(tmp_xy.first, tmp_xy.second))
						{
							mesh_path.push_back(tmp_xy);
							il_out << tmp_xy.first << " " << tmp_xy.second << endl;
							//prev_vid = curr_vid;
							//curr_vid = Round(tmp_xy.first)*vr_size_2 + Round(tmp_xy.second);
							curr = tmp_xy;
							cout << "This endPoint is a CriticalPoint！" << endl;
							break;
						}
						else
						{
							double x_dir = EvaluateGradX(tmp_xy.first, tmp_xy.second);
							double y_dir = EvaluateGradY(tmp_xy.first, tmp_xy.second);
							tmp_xy.first += x_dir / sqrt(norm(x_dir, y_dir));
							tmp_xy.second += y_dir / sqrt(norm(x_dir, y_dir));
							bool b = true;
							for (int t = mesh_path.size() - 4; t < mesh_path.size(); ++t)
							{
								if (tmp_xy.first - mesh_path[t].first>1 || tmp_xy.second - mesh_path[t].second > 1)
								{
									b = false;
								}
							}
							if (b == true)
								break;
						}
					}
					mesh_path.push_back(tmp_xy);
					il_out << tmp_xy.first << " " << tmp_xy.second << endl;
					//prev_vid = curr_vid;
					//curr_vid = Round(tmp_xy.first)*vr_size_2 + Round(tmp_xy.second);
					curr = tmp_xy;
					if (isMinimal(curr.first, curr.second))
					{
						cout << "ERROR,This line is ascend！" << endl;
						break;
					}
				}
				il.startIndex = it->xy_local;
				il.endIndex = mesh_path[mesh_path.size() - 1];
				il_out << endl;
			}
		}
		for (vector<CriticalPoint>::iterator it = msc.saddles.begin(); it != msc.saddles.end(); ++it)
		{
			if (it->type == SADDLE)
			{
				CriticalPoint& sad = *it;
				pair<double, double> curr(sad.xy_local);
				int prev_vid = -1, curr_vid = sad.meshIndex;
				msc.il_vec.push_back(IntegrationLine());
				IntegrationLine &il = msc.il_vec[msc.il_vec.size() - 1];
				PATH& mesh_path = il.path;
				mesh_path.push_back(sad.xy_local);
				il_out << sad.xy_local.first << " " << sad.xy_local.second << endl;
				while (!isMaximal(curr.first, curr.second))
				{
					if (isBoundary(curr.first, curr.second))
						break;
					while (isSaddle(curr.first, curr.second))
					{
						pair<double, double> saddle_positon = SaddlePosition(curr.first, curr.second);
						pair<double, double> xy = curr;
						double eig_vector_x = -msc.cp_vec[saddle_positon.first*vr_size_2 + saddle_positon.second].eig_vector2.first;
						double eig_vector_y = -msc.cp_vec[saddle_positon.first*vr_size_2 + saddle_positon.second].eig_vector2.second;
						xy = getTheSaddleBeginDirection(xy, make_pair(eig_vector_x, eig_vector_y));
						mesh_path.push_back(xy);
						il_out << xy.first << " " << xy.second << endl;
						//curr_vid = Round(xy.first)*vr_size_2 + Round(xy.second);
						curr = xy;
					}
					if (isBoundary(curr.first, curr.second))
						break;
					pair<double, double> tmp_xy;
					tmp_xy = getGradDirectionUp(curr);

					if (norm(tmp_xy.first - mesh_path[mesh_path.size() - 2].first, tmp_xy.second - mesh_path[mesh_path.size() - 2].second) < APPROACH*APPROACH || norm(tmp_xy.first - mesh_path[mesh_path.size() - 1].first, tmp_xy.second - mesh_path[mesh_path.size() - 1].second) < APPROACH*APPROACH)
					{
						if ((EvaluateGradX(tmp_xy.first, tmp_xy.second) < LARGE_ZERO_EPSILON&&
							EvaluateGradY(tmp_xy.first, tmp_xy.second) < LARGE_ZERO_EPSILON) ||
							isKeyPoint(tmp_xy.first, tmp_xy.second))
						{
							mesh_path.push_back(tmp_xy);
							il_out << tmp_xy.first << " " << tmp_xy.second << endl;
							//prev_vid = curr_vid;
							//curr_vid = Round(tmp_xy.first)*vr_size_2 + Round(tmp_xy.second);
							curr = tmp_xy;
							cout << "This endPoint is a CriticalPoint！" << endl;
							break;
						}
						else
						{
							double x_dir = EvaluateGradX(tmp_xy.first, tmp_xy.second);
							double y_dir = EvaluateGradY(tmp_xy.first, tmp_xy.second);
							tmp_xy.first += x_dir / sqrt(norm(x_dir, y_dir));
							tmp_xy.second += y_dir / sqrt(norm(x_dir, y_dir));
							bool b = true;
							for (int t = mesh_path.size() - 4; t < mesh_path.size(); ++t)
							{
								if (tmp_xy.first - mesh_path[t].first>1 || tmp_xy.second - mesh_path[t].second > 1)
								{
									b = false;
								}
							}
							if (b == true)
								break;
						}
					}
					mesh_path.push_back(tmp_xy);
					il_out << tmp_xy.first << " " << tmp_xy.second << endl;
					//prev_vid = curr_vid;
					//curr_vid = Round(tmp_xy.first)*vr_size_2 + Round(tmp_xy.second);
					curr = tmp_xy;
					if (isMinimal(curr.first, curr.second))
					{
						cout << "ERROR,This line is ascend！" << endl;
						break;
					}
				}
				il.startIndex = it->xy_local;
				il.endIndex = mesh_path[mesh_path.size() - 1];
				il_out << endl;
			}
		}
		for (vector<CriticalPoint>::iterator it = msc.saddles.begin(); it != msc.saddles.end(); ++it)
		{
			if (it->type == SADDLE)
			{
				CriticalPoint& sad = *it;
				pair<double, double> curr(sad.xy_local);
				int prev_vid = -1, curr_vid = sad.meshIndex;
				msc.il_vec.push_back(IntegrationLine());
				IntegrationLine &il = msc.il_vec[msc.il_vec.size() - 1];
				PATH& mesh_path = il.path;
				mesh_path.push_back(sad.xy_local);
				il_out << sad.xy_local.first << " " << sad.xy_local.second << endl;
				while (!isMaximal(curr.first, curr.second))
				{
					if (isBoundary(curr.first, curr.second))
						break;
					while (isSaddle(curr.first, curr.second))
					{
						pair<double, double> saddle_positon = SaddlePosition(curr.first, curr.second);
						pair<double, double> xy = curr;
						double eig_vector_x = msc.cp_vec[saddle_positon.first*vr_size_2 + saddle_positon.second].eig_vector2.first;
						double eig_vector_y = msc.cp_vec[saddle_positon.first*vr_size_2 + saddle_positon.second].eig_vector2.second;
						xy = getTheSaddleBeginDirection(xy, make_pair(eig_vector_x, eig_vector_y));
						mesh_path.push_back(xy);
						il_out << xy.first << " " << xy.second << endl;
						//curr_vid = Round(xy.first)*vr_size_2 + Round(xy.second);
						curr = xy;
					}
					if (isBoundary(curr.first, curr.second))
						break;
					pair<double, double> tmp_xy;
					tmp_xy = getGradDirectionUp(curr);

					if (norm(tmp_xy.first - mesh_path[mesh_path.size() - 2].first, tmp_xy.second - mesh_path[mesh_path.size() - 2].second) < APPROACH*APPROACH || norm(tmp_xy.first - mesh_path[mesh_path.size() - 1].first, tmp_xy.second - mesh_path[mesh_path.size() - 1].second) < APPROACH*APPROACH)
					{
						if ((EvaluateGradX(tmp_xy.first, tmp_xy.second) < LARGE_ZERO_EPSILON&&
							EvaluateGradY(tmp_xy.first, tmp_xy.second) < LARGE_ZERO_EPSILON) ||
							isKeyPoint(tmp_xy.first, tmp_xy.second))
						{
							mesh_path.push_back(tmp_xy);
							il_out << tmp_xy.first << " " << tmp_xy.second << endl;
							//prev_vid = curr_vid;
							//curr_vid = Round(tmp_xy.first)*vr_size_2 + Round(tmp_xy.second);
							curr = tmp_xy;
							cout << "This endPoint is a CriticalPoint！" << endl;
							break;
						}
						else
						{
							double x_dir = EvaluateGradX(tmp_xy.first, tmp_xy.second);
							double y_dir = EvaluateGradY(tmp_xy.first, tmp_xy.second);
							tmp_xy.first += x_dir / sqrt(norm(x_dir, y_dir));
							tmp_xy.second += y_dir / sqrt(norm(x_dir, y_dir));
							bool b = true;
							for (int t = mesh_path.size() - 4; t < mesh_path.size(); ++t)
							{
								if (tmp_xy.first - mesh_path[t].first>1 || tmp_xy.second - mesh_path[t].second > 1)
								{
									b = false;
								}
							}
							if (b == true)
								break;
						}
					}
					mesh_path.push_back(tmp_xy);
					il_out << tmp_xy.first << " " << tmp_xy.second << endl;
					//prev_vid = curr_vid;
					//curr_vid = Round(tmp_xy.first)*vr_size_2 + Round(tmp_xy.second);
					curr = tmp_xy;
					if (isMinimal(curr.first, curr.second))
					{
						cout << "ERROR,This line is ascend！" << endl;
						break;
					}
				}
				il.startIndex = it->xy_local;
				il.endIndex = mesh_path[mesh_path.size() - 1];
				il_out << endl;
			}
		}
		//ascend first part
		//for (vector<CriticalPoint>::iterator it = msc.saddles.begin(); it != msc.saddles.end(); ++it)
		//{			
		//		double px, py, posx, posy, dstepx, dstepy, testx, testy;								
		//		double pval, testval, step, l, dist,dir_x,dir_y;
		//		unsigned int path, i, j, num;
		//		
		//		for (path = 0; path < 2;++path)
		//		{
		//			msc.il_vec.push_back(IntegrationLine());
		//			IntegrationLine &il = msc.il_vec[msc.il_vec.size() - 1];
		//			PATH& mesh_path = il.path;
		//			num = 0; 
		//			posx = it->xy_local.first; posy = it->xy_local.second;
		//			mesh_path.push_back(make_pair(posx,posy));
		//			px = it->xy_local.first; py = it->xy_local.second;
		//			testx = it->xy_local.first; testy = it->xy_local.second;
		//			
		//			dir_x = it->eig_vector1.first;
		//			dir_y = it->eig_vector1.second;
		//			if (path>0)
		//			{
		//				dir_x = -dir_x;
		//				dir_y = -dir_y;
		//			}
		//			while (1)
		//			{
		//				posx = px + dir_x*MAXSTEP;
		//				posy = py + dir_y*MAXSTEP;
		//				if (posx<3||posx>vr_size_1-4||
		//					posy<3||posy>vr_size_2-4)
		//				{
		//					il.startIndex = it->xy_local;
		//					il.endIndex = mesh_path[mesh_path.size() - 1];
		//					break;
		//				}
		//				pval = EvaluateVal(px, py);
		//				step = MAXSTEP;
		//				while (step>=EPS)
		//				{
		//					posx = px+dir_x*step;
		//					posy = py + dir_y*step;
		//					testx = posx; testy = posy;
		//					testval = EvaluateVal(testx, testy);
		//					step *= 0.618;
		//					if (testval > pval)
		//					{
		//						posx = px + dir_x*step;
		//						posy = py + dir_y*step;
		//						px = posx; py = posy;
		//						if (pval >= EvaluateVal(px,py))
		//							px = testx; py = testy;
		//						dir_x = EvaluateGradX(px, py);
		//						dir_y = EvaluateGradY(px, py);
		//						//l = NORM2(dir_x, dir_y);
		//						//l = norm(dir_x, dir_y);
		//						//if (l > 0)
		//						//{
		//						//	//l = 1.0 /sqrt(l);									
		//						//	dir_x *= l;
		//						//	dir_y *= l;
		//						//}
		//					break;
		//					}
		//				}
		//				if (step < EPS)
		//				{
		//					dist = MAXSTEP;
		//					for (int i = 0; i < msc.maxPoints.size(); ++i)
		//					{
		//						dstepx = msc.maxPoints[i].xy_local.first - px;
		//						dstepy = msc.maxPoints[i].xy_local.second - py;
		//						if (norm(dstepx, dstepy) < dist)
		//						{
		//							j = i;
		//							dist = norm(dstepx,dstepy);
		//						}
		//					}
		//					if (dist<MAXSTEP)
		//					{
		//						posx = msc.maxPoints[j].xy_local.first;
		//						posy = msc.maxPoints[j].xy_local.second;
		//						mesh_path.push_back(make_pair(posx, posy));
		//						il.startIndex = it->xy_local;
		//						il.endIndex = mesh_path[mesh_path.size() - 1];
		//					}
		//					else
		//					{
		//						CriticalPoint newmax;
		//						newmax.xy_local = make_pair(posx, posy);
		//						newmax.type = CriticalPointType::MAXIMAL;
		//						msc.maxPoints.push_back(newmax);
		//						msc.keyPoint.push_back(newmax);
		//						mesh_path.push_back(make_pair(posx, posy));
		//						il.startIndex = it->xy_local;
		//						il.endIndex = mesh_path[mesh_path.size() - 1];
		//					}
		//					break;
		//				}
		//				else
		//				{
		//					dstepx = mesh_path[num].first-posx;
		//					dstepy = mesh_path[num].second - py;
		//					if (norm(dstepx, dstepy) > MAXSTEP*MAXSTEP)
		//					{
		//						mesh_path.push_back(make_pair(posx, posy));
		//						++num;
		//					}
		//				}
		//			}
		//		}
		//	}		
  		return true;
	}
	bool ILTracer::isMaximal(double x, double y)
	{
		int index = Round(x)*vr_size_2 + Round(y);
		for (int i = -1; i < 2; ++i)
		{
			for (int j = -1; j < 2; ++j)
			{
				if (msc.cp_vec[index - j*vr_size_2 - i].type == MAXIMAL)
					return true;
			}
		}
		return 	false;
	}
	bool ILTracer::isMinimal(double x, double y)
	{
		int index = Round(x)*vr_size_2 + Round(y);
		for (int i = -1; i < 2; ++i)
		{
			for (int j = -1; j < 2; ++j)
			{
				if (msc.cp_vec[index - j*vr_size_2 - i].type == MINIMAL)
					return true;
			}
		}
		return 	false;
	}
	bool ILTracer::isSaddle(double x, double y)
	{
		int index = Round(x)*vr_size_2 + Round(y);
		for (int i = -1; i < 2; ++i)
		{
			for (int j = -1; j < 2; ++j)
			{
				if (msc.cp_vec[index - j*vr_size_2 - i].type == SADDLE)
					return true;
			}
		}
		return 	false;		
	}
	bool ILTracer::isKeyPoint(double x, double y)
	{
		if (isMaximal(x, y) || isMinimal(x, y) || isSaddle(x, y))
			return true;
		else
			return false;
	}
	pair<double,double> ILTracer::SaddlePosition(double x, double y)
	{
		pair<double, double> sad_point;
		int index = Round(x)*vr_size_2 + Round(y);
		for (int i = -1; i < 2; ++i)
		{
			for (int j = -1; j < 2; ++j)
			{
				if (msc.cp_vec[index - j*vr_size_2 - i].type == SADDLE)
					sad_point = msc.cp_vec[index - j*vr_size_2 - i].xy_local;
			}
		}
		return sad_point;
	}
	bool ILTracer::isBoundary(double x, double y)
	{
		if (x >= vr_size_1 - 4 || y >= vr_size_2 - 4 || x < 3 || y < 3)
			return true;
		return false;
	}
	int ILTracer::Round(double r)
	{
		
		//return r > 0 ? ceil(r) : floor(r);
		return floor(r+0.5);
	}
	pair<double, double>ILTracer::getTheSaddleBeginDirection(pair<double, double> xy, pair<double, double> eig_vector)
	{
		double cur_x = xy.first;
		double cur_y = xy.second;
		double eig_vector_x = eig_vector.first;
		double eig_vector_y = eig_vector.second;

		return make_pair(cur_x + eig_vector_x, cur_y + eig_vector_y);
	}
	pair<double, double> ILTracer::getGradDirectionDown(pair<double, double> xy)
	{
		//四阶龙格库塔法 ,待检验，Round向下取整

		double x = xy.first, y = xy.second;
		double k1x = EvaluateGradX(x, y);
		//double k1x = msc.cp_vec[x*vr_size_2 + y].dif.first;
		double k2x = EvaluateGradX(x + k1x/2, y + 1/2);
		double k3x = EvaluateGradX(x + k2x/2, y + 1/2);
		double k4x = EvaluateGradX(x + k3x, y + 1);
		double tmp_x = (k1x + 2 * k2x + 2 * k3x + k4x) / 6;

		//double k1y = msc.cp_vec[x*vr_size_2 + y].dif.second;
		double k1y = EvaluateGradY(x, y);
		double k2y = EvaluateGradY(x + 1/2, y + k1y/2);
		double k3y = EvaluateGradY(x + 1/2, y + k2y/2);
		double k4y = EvaluateGradY(x + 1, y + k3y);
		double tmp_y = (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
		/*int k2x = Round(msc.cp_vec[(x + k1x)*vr_size_2 + y + 1].dif.first);
		int k3x = Round(msc.cp_vec[(x + k2x)*vr_size_2 + y + 1].dif.first);
		int k1y = Round(msc.cp_vec[x*vr_size_2 + y].dif.second);
		int k2y = Round(msc.cp_vec[(x + 1)*vr_size_2 + y + k1y].dif.second);
		int k3y = Round(msc.cp_vec[(x + 1)*vr_size_2 + y + k2y].dif.second);
		int tmp_x = Round((msc.cp_vec[x*vr_size_2 + y].dif.first +
		2 * msc.cp_vec[(x + k1x)*vr_size_2 + y + 1].dif.first +
		2 * msc.cp_vec[(x + k2x)*vr_size_2 + y + 1].dif.first +
		msc.cp_vec[(x + k3x)*vr_size_2 + y + 2].dif.first) / 3);
		int tmp_y = Round((msc.cp_vec[x*vr_size_2 + y].dif.second +
		2 * msc.cp_vec[(x + 1)*vr_size_2 + y + k1y].dif.second +
		2 * msc.cp_vec[(x + 1)*vr_size_2 + y + k2y].dif.second +
		msc.cp_vec[(x + 2)*vr_size_2 + y + k3y].dif.second) / 3);*/
		
		double next_X = x - tmp_x / sqrt(norm(tmp_x, tmp_y));
		double next_Y = y - tmp_y/ sqrt(norm(tmp_x, tmp_y));;

		return make_pair(next_X, next_Y);
	}
	pair<double, double> ILTracer::getGradDirectionUp(pair<double, double> xy)
	{
		//四阶龙格库塔法 ,待检验，Round向下取整

		double x = xy.first, y = xy.second;
		double k1x = EvaluateGradX(x, y);

		double k2x = EvaluateGradX(x + k1x / 2, y + 1 / 2);
		double k3x = EvaluateGradX(x + k2x / 2, y + 1 / 2);
		double k4x = EvaluateGradX(x + k3x, y + 1);
		double tmp_x = (k1x + 2 * k2x + 2 * k3x + k4x) / 6;

		double k1y = EvaluateGradY(x, y);
		double k2y = EvaluateGradY(x + 1 / 2, y + k1y / 2);
		double k3y = EvaluateGradY(x + 1 / 2, y + k2y / 2);
		double k4y = EvaluateGradY(x + 1, y + k3y);
		double tmp_y = (k1y + 2 * k2y + 2 * k3y + k4y) / 6;

		double next_X = x + tmp_x / sqrt(norm(tmp_x, tmp_y));
		double next_Y = y + tmp_y / sqrt(norm(tmp_x, tmp_y));;

		return make_pair(next_X, next_Y);
	}
	
}

