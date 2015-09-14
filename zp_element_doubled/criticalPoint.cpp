#include <iostream>
#include "zp_element_doubled.h"
#include "box_spline.h"
#include "point_evaluator_6direction2.h"
#include "point_evaluator_7direction.h"
#include "mscomplex.h"
#include "criticalPoint.h"
#include "vtkReader.h"
#include <map>
namespace msc2d
{											//2ddata  gd1wie220_20  large=0.1  5 0.01
	const double LARGE_ZERO_EPSILON = 0.04;//401数据阈值0.0025,2e-6;taiwan_2_double3D_10_filter_50:5e-4,1e-100
	const double SMALL_ZERO_EPSILON = 1e-100;
	using namespace std;
	vtkReader vr;
	point_evaluator_7direction* p7;

	CriticalPointFind::CriticalPointFind(MSComplex2D& _msc) :msc(_msc){}
	CriticalPointFind::~CriticalPointFind(){}
	int cut(int i, int N)
	{
		int tmp = ((i<0) ? 0 : i);
		return ((tmp>(N - 1)) ? (N - 1) : tmp);
	}
	void CriticalPointFind::findCriticalPoint()
	{
		cout << "Find Critical Point Begin!" << endl;
		box_spline *bs;
		bs = new box_spline();
		vr.loadFile("D:\\newData\\topo\\gd1wie220_20.nak");//201数据：由于201数据步长为0.1，求偏导数的默认步长为1，故求出结果为实际结果的十分之一；
		double x_ordinates[10000], y_ordinates[10000];

		int m = 0;
		int position[2];
		double sum, dif__x, dif__y, dif__yy, dif__xx, dif__xy, Value, dif__x1, dif__y1;
		double distance[2];
		double diff_x, diff_y;
		double f_with_diff, f_with_diff2;
		double value1, value3, value2, value4, value5, value0;
		for (int i = 0; i < vr.dim[0]; i++)
			x_ordinates[i] = i;
		for (int i = 0; i < vr.dim[1]; i++)
			y_ordinates[i] = i;
		CriticalPointArray &cp_vec = msc.cp_vec;//一维存放,cp_vec设为public类型可以访问，private不能访问,当定义为friend类时，可以访问私有类型
		CriticalPointArray &minPoint = msc.minPoints;
		CriticalPointArray &maxPoint = msc.maxPoints;
		CriticalPointArray &saddles = msc.saddles;
		CriticalPointArray &keyPoint = msc.keyPoint;

		cp_vec.clear(); minPoint.clear(); maxPoint.clear(); saddles.clear();
		
		size_t k = 0;
		double gxmin, gymin,  gxmax, gymax, tgx, tgy;
		for (int i = 0; i < vr.dim[0]; i++)
		{
			for (int j = 0; j < vr.dim[1]; j++)
			{
				sum = 0; dif__x = 0; dif__y = 0; dif__yy = 0, dif__xy = 0, Value = 0, dif__xx = 0,dif__x1=0,dif__y1=0;
				CriticalPoint cp;
				cp.xy_local.first = i;
				cp.xy_local.second = j;
				cp.meshIndex = i*vr.dim[1] + j;				

				gxmin = gymin = 1; gxmax = gymax = -1;
				for (int m1 = 0; m1 < 7;++m1)
				{
					for (int m2 = 0; m2 < 7;++m2)
					{
						tgx = vr.getData(cut(i + m2 + 1, vr.dim[0]), cut(j + m1, vr.dim[1])) - vr.getData(cut(i + m2 , vr.dim[0]), cut(j + m1, vr.dim[1]));
						tgy = vr.getData(cut(i + m2 , vr.dim[0]), cut(j + m1+1, vr.dim[1])) - vr.getData(cut(i + m2, vr.dim[0]), cut(j + m1, vr.dim[1]));
						if (tgx < gxmin)
							gxmin = tgx;
						if (tgy < gymin)
							gymin = tgy;
						if (tgx>gxmax)
							gxmax = tgx;
						if (tgy>gymax)
							gymax = tgy;
					}
				}
				if (gxmin*gxmax < 0 && gymax*gymin < 0)
				{
					bool flag = true;
					//求导
					if (i<3 || j<3 || i>vr.dim[0] - 4 || j>vr.dim[1] - 4)
					{
						dif__x = vr.getData(cut(floor(i + 0.5) + 1, vr.dim[0]), j) - vr.getData(floor(i + 0.5), floor(j + 0.5));
						dif__y = vr.getData(floor(i + 0.5), cut(floor(j + 0.5) + 1, vr.dim[1])) - vr.getData(floor(i + 0.5), floor(j + 0.5));
					}
					else
					{
						for (int p = 0; p < 7; p++)
						{
							for (int q = 0; q < 6; q++)
							{
								position[0] = cut((int)floor(x_ordinates[i] + 0.5) + p - 3, vr.dim[0]);
								position[1] = cut((int)floor(y_ordinates[j] + 0.5) + q - 3, vr.dim[1]);
								distance[0] = position[0] - x_ordinates[i];
								distance[1] = position[1] - y_ordinates[j];

								value1 = bs->compute_gradient_x(distance);
								diff_x = vr.getData(cut(position[0] + 1, vr.dim[0]), position[1])
									+ vr.getData(cut(position[0] - 1, vr.dim[0]), position[1])
									- 2 * vr.getData(position[0], position[1]);
								diff_y = vr.getData(position[0], cut(position[1] + 1, vr.dim[1]))
									+ vr.getData(position[0], cut(position[1] - 1, vr.dim[1]))
									- 2 * vr.getData(position[0], position[1]);

								f_with_diff = vr.getData(position[0], position[1]) - 1.0*(diff_x + diff_y) / 8;
								dif__x += f_with_diff*value1;
								//dif__x += vr.getData(position[0], position[1])*value1;
							}
						}
						for (int p = 0; p < 6; p++)
						{
							for (int q = 0; q < 7; q++)
							{
								position[0] = cut((int)floor(x_ordinates[i] + 0.5) + p - 3, vr.dim[0]);
								position[1] = cut((int)floor(y_ordinates[j] + 0.5) + q - 3, vr.dim[1]);
								distance[0] = position[0] - x_ordinates[i];
								distance[1] = position[1] - y_ordinates[j];

								value2 = bs->compute_gradient_y(distance);
								diff_x = vr.getData(cut(position[0] + 1, vr.dim[0]), position[1])
									+ vr.getData(cut(position[0] - 1, vr.dim[0]), position[1])
									- 2 * vr.getData(position[0], position[1]);
								diff_y = vr.getData(position[0], cut(position[1] + 1, vr.dim[1]))
									+ vr.getData(position[0], cut(position[1] - 1, vr.dim[1]))
									- 2 * vr.getData(position[0], position[1]);

								f_with_diff = vr.getData(position[0], position[1]) - 1.0*(diff_x + diff_y) / 8;
								dif__y += f_with_diff*value2;
								//dif__y += vr.getData(position[0], position[1])*value2;

							}
						}
						for (int w1 = i - 3; w1 < i; ++w1)
						{
							for (int w2 = j - 3; w2 < j; ++w2)
							{
								if (cp_vec[w1*vr.dim[1] + w2].type != REGULAR)
								{
									flag = false;
								}
							}
						}
					}

					//dif__x1 = EvaluateGradX(i, j);					
					//dif__y1 = EvaluateGradY(i, j);
					if (dif__x< LARGE_ZERO_EPSILON &&dif__x>-LARGE_ZERO_EPSILON)
						dif__x = 0;
					if (dif__y<LARGE_ZERO_EPSILON&&dif__y>-LARGE_ZERO_EPSILON)
						dif__y = 0;
					cp.dif = make_pair(dif__x, dif__y);
					if (flag&&dif__x == 0 && dif__y == 0&&i>2&&j>2)
					{						
						double eig_value[2], eig_vector[2][2], dif[2][2], eps = 1e-50;
						const int Dim = 2, nJt = 1e20;
						dif[0][0] = EvaluateGradXX(i, j); 
						dif[0][1] = EvaluateGradXY(i, j);
						dif[1][0] = EvaluateGradXY(i, j);
						dif[1][1] = EvaluateGradYY(i, j);
						//cout << i << " " << j << endl;
						if (JacbiCor(*dif, *eig_vector, eig_value, eps, nJt))
						{
							if (eig_value[1]<0 && eig_value[0]<0)
							{
								cp.type = MAXIMAL;
								cp.eig_vector1 = make_pair(0, 0);
								cp.eig_vector2 = make_pair(0, 0);						
								maxPoint.push_back(cp);
								//Value = 0;
							}
							else if (eig_value[0]>0 && eig_value[1]>0)
							{
								cp.type = MINIMAL;
								cp.eig_vector1 = make_pair(0, 0);
								cp.eig_vector2 = make_pair(0, 0);
								minPoint.push_back(cp);
								//Value = 125;
							}
							else
							{
								cp.type = SADDLE;
								//特征向量
								cp.eig_vector1 = make_pair(eig_vector[0][0], eig_vector[1][0]);
								cp.eig_vector2 = make_pair(eig_vector[0][1], eig_vector[1][1]);
								saddles.push_back(cp);
								//Value = 255;
							}
							//value <<i << " "<<j<<" "<<cp.type<<endl;
						}
						keyPoint.push_back(cp);
					}
					else
						cp.type = REGULAR;
					
				}	
				else
				{
					cp.type = REGULAR;
				}
				cp_vec.push_back(cp);														
			}		
		}
		delete bs;
		cout << "Critical point number:" << keyPoint.size() << endl;
		cout << "Max point number:" << maxPoint.size() << endl;
		cout << "Saddle point number:" << saddles.size() << endl;
		cout << "Min point number:" << minPoint.size() << endl;

 		cout << "End !" << endl << endl;
	}
	double EvaluateVal(const double x, const double y)
	{
		box_spline* bs;
		bs = new box_spline();
		int position[2];
		double tmpValue, val = 0;
		double distance[2];

		if (x<3 || y<3 || x>vr.dim[0] - 4 || y>vr.dim[1] - 4)
		{
			return vr.getData(floor(x + 0.5), floor(y + 0.5));
		}
		else
		{
			for (int p = 0; p < 7; p++)
			{
				for (int q = 0; q < 7; q++)
				{
					position[0] = cut((int)floor(x + 0.5) + p - 3, vr.dim[0]);
					position[1] = cut((int)floor(y + 0.5) + q - 3, vr.dim[1]);
					distance[0] = position[0] - x;
					distance[1] = position[1] - y;
					tmpValue = bs->compute_value(distance);
					val+=vr.getData(position[0], position[1])*tmpValue;
				}
			}
		}
		delete bs;

		return val;
	}
	double EvaluateGradX(const double x, const double y)
	{
		box_spline* bs;
		bs = new box_spline();
		int position[2];
		double tmpValue, val = 0, diff_x=0, diff_y=0, f_with_diff=0;
		double distance[2];

		if (x<3 || y<3 || x>vr.dim[0] - 4 || y>vr.dim[1] - 4)
		{
			return vr.getData(cut(floor(x + 0.5) + 1, vr.dim[0]),floor(y + 0.5)) - vr.getData(floor(x + 0.5), floor(y + 0.5));
		}
		else
		{
			for (int p = 0; p < 7; p++)
			{
				for (int q = 0; q < 6; q++)
				{
					position[0] = cut((int)floor(x + 0.5) + p - 3, vr.dim[0]);
					position[1] = cut((int)floor(y + 0.5) + q - 3, vr.dim[1]);
					distance[0] = position[0] - x;
					distance[1] = position[1] - y;

					tmpValue = bs->compute_gradient_x(distance);
					diff_x = vr.getData(cut(position[0] + 1, vr.dim[0]), position[1])
						+ vr.getData(cut(position[0] - 1, vr.dim[0]), position[1])
						- 2 * vr.getData(position[0], position[1]);
					diff_y = vr.getData(position[0], cut(position[1] + 1, vr.dim[1]))
						+ vr.getData(position[0], cut(position[1] - 1, vr.dim[1]))
						- 2 * vr.getData(position[0], position[1]);

					f_with_diff = vr.getData(position[0], position[1]) - 1.0*(diff_x + diff_y) / 8;
					val += f_with_diff*tmpValue;
				}
			}
		}
		delete bs;
		return val;
	}
	double EvaluateGradX1(const double x, const double y)
	{
		box_spline* bs;
		bs = new box_spline();
		int position[2];
		double tmpValue, val = 0, diff_x = 0, diff_y = 0, f_with_diff = 0;
		double distance[2];

		if (x<3 || y<3 || x>vr.dim[0] - 4 || y>vr.dim[1] - 4)
		{
			return vr.getData(cut(floor(x + 0.5) + 1, vr.dim[0]), floor(y + 0.5)) - vr.getData(floor(x + 0.5), floor(y + 0.5));
		}
		else
		{
			for (int p = 0; p < 7; p++)
			{
				for (int q = 0; q < 6; q++)
				{
					position[0] = cut((int)floor(x + 0.5) + p - 3, vr.dim[0]);
					position[1] = cut((int)floor(y + 0.5) + q - 3, vr.dim[1]);
					distance[0] = position[0] - x;
					distance[1] = position[1] - y;

					tmpValue = bs->compute_gradient_x(distance);
					val += vr.getData(position[0], position[1])*tmpValue;
				}
			}
		}
		delete bs;
		return val;
	}
	double EvaluateGradY(const double x, const double y)
	{
		box_spline* bs;
		bs = new box_spline();
		int position[2];
		double tmpValue, val = 0, diff_x = 0, diff_y = 0, f_with_diff = 0;
		double distance[2];

		if (x<3 || y<3 || x>vr.dim[0] - 4 || y>vr.dim[1] - 4)
		{
			return vr.getData(floor(x + 0.5), cut(floor(y + 0.5)+1,vr.dim[1])) - vr.getData(floor(x + 0.5), floor(y + 0.5));
		}
		else
		{
			for (int p = 0; p < 6; p++)
			{
				for (int q = 0; q < 7; q++)
				{												
					position[0] = cut((int)floor(x + 0.5) + p - 3, vr.dim[0]);
					position[1] = cut((int)floor(y + 0.5) + q - 3, vr.dim[1]);
					distance[0] = position[0] - x;
					distance[1] = position[1] - y;
					tmpValue = bs->compute_gradient_y(distance);
					diff_x = vr.getData(cut(position[0] + 1, vr.dim[0]), position[1])
						+ vr.getData(cut(position[0] - 1, vr.dim[0]), position[1])
						- 2 * vr.getData(position[0], position[1]);
					diff_y = vr.getData(position[0], cut(position[1] + 1, vr.dim[1]))
						+ vr.getData(position[0], cut(position[1] - 1, vr.dim[1]))
						- 2 * vr.getData(position[0], position[1]);

					f_with_diff = vr.getData(position[0], position[1]) - 1.0*(diff_x + diff_y) / 8;
					val += f_with_diff*tmpValue;
				}
			}
		}
		delete bs;
		return val;
	}
	double EvaluateGradXX(const double x, const double y)
	{
		box_spline* bs;
		bs = new box_spline();
		int position[2];
		double tmpValue, val = 0;
		double distance[2];

		if (x<3 || y<3 || x>vr.dim[0] - 4 || y>vr.dim[1] - 4)
		{
			return EvaluateGradX(cut(x + 1, vr.dim[0]), y) - EvaluateGradX(x, y);
		}
		else
		{
			for (int p = 0; p < 7; p++)
			{
				for (int q = 0; q < 5; q++)
				{
					position[0] = cut((int)floor(x + 0.5) + p - 3, vr.dim[0]);
					position[1] = cut((int)floor(y + 0.5) + q - 3, vr.dim[1]);
					distance[0] = position[0] - x;
					distance[1] = position[1] - y;
					tmpValue = bs->compute_gradient_xx(distance);
					val += vr.getData(position[0], position[1])*tmpValue;
				}
			}
		}
		delete bs;

		return val;
	}
	double EvaluateGradYY(const double x, const double y)
	{
		box_spline* bs;
		bs = new box_spline();
		int position[2];
		double tmpValue, val = 0;
		double distance[2];

		if (x<3 || y<3 || x>vr.dim[0] - 4 || y>vr.dim[1] - 4)
		{
			return EvaluateGradY(x, cut(y + 1, vr.dim[1])) - EvaluateGradY(x, y);
		}
		else
		{
			for (int p = 0; p < 5; p++)
			{
				for (int q = 0; q < 7; q++)
				{
					position[0] = cut((int)floor(x + 0.5) + p - 3, vr.dim[0]);
					position[1] = cut((int)floor(y + 0.5) + q - 3, vr.dim[1]);
					distance[0] = position[0] - x;
					distance[1] = position[1] - y;
					tmpValue = bs->compute_gradient_yy(distance);
					val += vr.getData(position[0], position[1])*tmpValue;
				}
			}
		}
		delete bs;

		return val;
	}
	double EvaluateGradXY(const double x, const double y)
	{
		box_spline* bs;
		bs = new box_spline();
		int position[2];
		double tmpValue, val = 0;
		double distance[2];

		if (x<3 || y<3 || x>vr.dim[0] - 4 || y>vr.dim[1] - 4)
		{
			return (EvaluateGradYY(x, y) + EvaluateGradXX(x, y))/2;
		}
		else
		{
			for (int p = 0; p < 6; p++)
			{
				for (int q = 0; q < 6; q++)
				{
					position[0] = cut((int)floor(x + 0.5) + p - 3, vr.dim[0]);
					position[1] = cut((int)floor(y + 0.5) + q - 3, vr.dim[1]);
					distance[0] = position[0] - x;
					distance[1] = position[1] - y;
					tmpValue = bs->compute_gradient_xy(distance);
					val += vr.getData(position[0], position[1])*tmpValue;
				}
			}
		}
		delete bs;

		return val;
	}
	bool CriticalPointFind::JacbiCor(double * pMatrix, double *pdblVects, double *pdbEigenValues, double dbEps, int nJt)
	{
		if (pMatrix[0] == pMatrix[3])
		{
			double eig1[2] = { 1, 0 }, eig2[2] = { 0, 1 };

			if (pMatrix[1] == 0 && pMatrix[2] == 0)
			{
				pdbEigenValues[0] = pMatrix[0];
				pdbEigenValues[1] = pMatrix[3];
				pdblVects[0] = *eig1;
				pdblVects[1] = *eig2;
				return true;
			}
		}
		const int nDim = 2;
		for (int i = 0; i < nDim; i++)
		{
			pdblVects[i*nDim + i] = 1.0f;
			for (int j = 0; j < nDim; j++)
			{
				if (i != j)
					pdblVects[i*nDim + j] = 0.0f;
			}
		}

		int nCount = 0;		//迭代次数
		while (1)
		{
			//在pMatrix的非对角线上找到最大元素
			double dbMax = pMatrix[1];
			int nRow = 0;
			int nCol = 1;
			for (int i = 0; i < nDim; i++)			//行
			{
				for (int j = 0; j < nDim; j++)		//列
				{
					double d = fabs(pMatrix[i*nDim + j]);

					if ((i != j) && (d> dbMax))
					{
						dbMax = d;
						nRow = i;
						nCol = j;
					}
				}
			}

			if (dbMax < dbEps)     //精度符合要求 
				break;

			if (nCount > nJt)       //迭代次数超过限制
				break;

			nCount++;

			double dbApp = pMatrix[nRow*nDim + nRow];
			double dbApq = pMatrix[nRow*nDim + nCol];
			double dbAqq = pMatrix[nCol*nDim + nCol];

			//计算旋转角度
			double dbAngle = 0.5*atan2(-2 * dbApq, dbAqq - dbApp);
			double dbSinTheta = sin(dbAngle);
			double dbCosTheta = cos(dbAngle);
			double dbSin2Theta = sin(2 * dbAngle);
			double dbCos2Theta = cos(2 * dbAngle);

			pMatrix[nRow*nDim + nRow] = dbApp*dbCosTheta*dbCosTheta +
				dbAqq*dbSinTheta*dbSinTheta + 2 * dbApq*dbCosTheta*dbSinTheta;
			pMatrix[nCol*nDim + nCol] = dbApp*dbSinTheta*dbSinTheta +
				dbAqq*dbCosTheta*dbCosTheta - 2 * dbApq*dbCosTheta*dbSinTheta;
			pMatrix[nRow*nDim + nCol] = 0.5*(dbAqq - dbApp)*dbSin2Theta + dbApq*dbCos2Theta;
			pMatrix[nCol*nDim + nRow] = pMatrix[nRow*nDim + nCol];

			for (int i = 0; i < nDim; i++)
			{
				if ((i != nCol) && (i != nRow))
				{
					int u = i*nDim + nRow;	//p  
					int w = i*nDim + nCol;	//q
					dbMax = pMatrix[u];
					pMatrix[u] = pMatrix[w] * dbSinTheta + dbMax*dbCosTheta;
					pMatrix[w] = pMatrix[w] * dbCosTheta - dbMax*dbSinTheta;
				}
			}

			for (int j = 0; j < nDim; j++)
			{
				if ((j != nCol) && (j != nRow))
				{
					int u = nRow*nDim + j;	//p
					int w = nCol*nDim + j;	//q
					dbMax = pMatrix[u];
					pMatrix[u] = pMatrix[w] * dbSinTheta + dbMax*dbCosTheta;
					pMatrix[w] = pMatrix[w] * dbCosTheta - dbMax*dbSinTheta;
				}
			}

			//计算特征向量
			for (int i = 0; i < nDim; i++)
			{
				int u = i*nDim + nRow;		//p   
				int w = i*nDim + nCol;		//q
				dbMax = pdblVects[u];
				pdblVects[u] = pdblVects[w] * dbSinTheta + dbMax*dbCosTheta;
				pdblVects[w] = pdblVects[w] * dbCosTheta - dbMax*dbSinTheta;
			}

		}

		//对特征值进行排序以及重新排列特征向量,特征值即pMatrix主对角线上的元素
		std::map<double, int> mapEigen;
		for (int i = 0; i < nDim; i++)
		{
			pdbEigenValues[i] = pMatrix[i*nDim + i];

			mapEigen.insert(make_pair(pdbEigenValues[i], i));
		}

		double *pdbTmpVec = new double[nDim*nDim];
		std::map<double, int>::reverse_iterator iter = mapEigen.rbegin();
		for (int j = 0; iter != mapEigen.rend(), j < nDim; ++iter, ++j)
		{
			for (int i = 0; i < nDim; i++)
			{
				pdbTmpVec[i*nDim + j] = pdblVects[i*nDim + iter->second];
			}

			//特征值重新排列
			pdbEigenValues[j] = iter->first;
		}

		//设定正负号
		for (int i = 0; i < nDim; i++)
		{
			double dSumVec = 0;
			for (int j = 0; j < nDim; j++)
				dSumVec += pdbTmpVec[j * nDim + i];
			if (dSumVec < 0)
			{
				for (int j = 0; j < nDim; j++)
					pdbTmpVec[j * nDim + i] *= -1;
			}
		}

		memcpy(pdblVects, pdbTmpVec, sizeof(double)*nDim*nDim);
		delete[]pdbTmpVec;

		return 1;
	} 
}
