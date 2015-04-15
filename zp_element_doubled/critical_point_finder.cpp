#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "vtkReader.h"
#include "zp_element_doubled.h"
#include "box_spline.h"
#include "point_evaluator_6direction2.h"
#include "point_evaluator_7direction.h"
#include <utility>
#include "critical_point_finder.h"
#include "mscomplex.h"
#include <map>
namespace msc2d
{
	using namespace std;
	const double LARGE_ZERO_EPSILON = 1e-14;//401sin 阈值0.0025,2e-6;
	const double SMALL_ZERO_EPSILON = 1e-50;
	using namespace std;
	//vtkReader vr;
	
	 VTK vtk;
	point_evaluator_7direction* p7;
	
	void reconstructor::set_dimension(int L, int W)
	{
		this->L = L;
		this->W = W;
	}

	void reconstructor::set_step(double stepx, double stepy)
	{
		this->stepx = stepx;
		this->stepy = stepy;

		this->nx = (int)floor(((this->L) - 1) / this->stepx) + 1;
		this->ny = (int)floor(((this->W) - 1) / this->stepy) + 1;
	}

	int cut(int i, int N)
	{
		int tmp = ((i<0) ? 0 : i);
		return ((tmp>(N - 1)) ? (N - 1) : tmp);
	}

	CPFinder::CPFinder(MSComplex2D& _msc) :msc(_msc){}
	CPFinder::~CPFinder(){}
	

	void CPFinder::findCriticalPoints()
	{
		
		string path = "D:/401data/sin401.vtk";
		vtk.Load(path, 0);
		int dimX = vtk.dimX, dimY = vtk.dimY;
		double *data = vtk.data;
		
		vtk.GaussianFilter(3, 0.5, 0.01,data,dimX,dimY);
		

		box_spline* bs;
		bs = new box_spline();
		//FILE * fp;
		
		//vr.loadFile("D:\\topology2ddata\\topology2d.vtk");//201数据：由于201数据步长为0.1，求偏导数的默认步长为1，故求出结果为实际结果的十分之一；
		double x_ordinates[10000], y_ordinates[10000];//第一个参数为列，按列读取数据
		ofstream value("D:\\401data\\sin401_cp_filter.txt");
		/*ofstream dif_x("D:\\topology2ddata\\topology2d_dif_x.txt");
		ofstream dif_y("D:\\topology2ddata\\topology2d_dif_y.txt");
		ofstream dif_xx("D:\\topology2ddata\\topology2d_dif_xx.txt");
		ofstream dif_xy("D:\\topology2ddata\\topology2d_dif_xy.txt");
		ofstream dif_yy("D:\\topology2ddata\\topology2d_dif_yy.txt");*/
		//cout<<vr.getData(0,0)<<" "<<vr.getData(0,1)<<" "<<vr.getData(0,2)<<endl;
		//cout<<vr.getData(1,0)<<" "<<vr.getData(1,1)<<" "<<vr.getData(1,2)<<endl;

		int m = 0;
		int position[2];
		double sum, dif__x, dif__y, dif__yy, dif__xx, dif__xy, Value;
		double distance[2];
		double diff_x, diff_y;
		double f_with_diff, f_with_diff2;
		double value1, value3, value2, value4, value5, value0;
		for (int i = 0; i < dimX; i++)
			x_ordinates[i] = i;
		for (int i = 0; i < dimY; i++)
			y_ordinates[i] = i;
		
		cout << "Find critical point begin!" << endl;
		CriticalPointArray &cp_vec = msc.cp_vec;//一维存放,cp_vec设为public类型可以访问，private不能访问,当定义为friend类时，可以访问私有类型
		CriticalPointArray &minPoint = msc.minPoints;
		CriticalPointArray &maxPoint = msc.maxPoints;
		CriticalPointArray &saddles = msc.saddles;
		CriticalPointArray &keyPoint = msc.keyPoint;

		cp_vec.clear(); minPoint.clear(); maxPoint.clear(); saddles.clear();
		//cp_vec.resize(dimX * dimY);
		size_t k = 0;

		for (int i = 0; i < dimY; i++)
		{
			for (int j = 0; j < dimX; j++)
			{
				sum = 0; dif__x = 0; dif__y = 0; dif__yy = 0, dif__xy = 0, Value = 0, dif__xx = 0;
				//value
				
				if (i <3 ||  j <3  || i >=dimX - 3 ||  j >= dimY - 3 )		
					Value = data[i*dimX + j];//Value = vr.getData(i, j);					
				else
				{
					for (int p = 0; p < 7; p++)
					{
						for (int q = 0; q < 7; q++)
						{
							position[0] = cut((int)floor(x_ordinates[j] + 0.5) + p - 3, dimX);
							position[1] = cut((int)floor(y_ordinates[i] + 0.5) + q - 3, dimY);
							distance[0] = position[0] - x_ordinates[j];
							distance[1] = position[1] - y_ordinates[i];
							value0 = bs->compute_value(distance);
							diff_x = data[cut(position[0] + 1, dimX)*dimX + position[1]]
								+ data[cut(position[0] - 1, dimX)*dimX + position[1]]
								- 2 * data[position[0] * dimX + position[1]  ];
							diff_x = data[position[0] * dimX + cut(position[1] + 1, dimY)]
								+ data[position[0] * dimX + cut(position[1] - 1, dimY) ]
								- 2 * data[position[0] * dimX + position[1]  ];
							
							f_with_diff = data[position[0] * dimX + position[1]] - 1.0*(diff_x + diff_y) / 8;
							Value += f_with_diff*value0;
							
						}
					}
				}
				if (Value > 255)Value = 255;
				if (Value < 0)Value = 0;
				

				//求导
				for (int p = 0; p < 7; p++)
				{
					for (int q = 0; q < 6; q++)
					{
						position[0] = cut((int)floor(x_ordinates[j] + 0.5) + p - 3, dimX);
						position[1] = cut((int)floor(y_ordinates[i] + 0.5) + q - 3, dimY);
						distance[0] = position[0] - x_ordinates[j]+0.5;
						distance[1] = position[1] - y_ordinates[i];

						value1 = bs->compute_gradient_x(distance);								
						/*diff_x = vr.getData(cut(position[0] + 1, dimX), position[1])
							+ vr.getData(cut(position[0] - 1, dimX), position[1])
							- 2 * vr.getData(position[0], position[1]);
						diff_y = vr.getData(position[0], cut(position[1] + 1, dimY))
							+ vr.getData(position[0], cut(position[1] - 1, dimY))
							- 2 * vr.getData(position[0], position[1]);

						f_with_diff = vr.getData(position[0], position[1]) - 1.0*(diff_x + diff_y) / 8;
						dif__x += f_with_diff*value1;*/
						dif__x += data[position[0] * dimX + position[1]]*value1;
					}
				}
				for (int p = 0; p < 6; p++)
				{
					for (int q = 0; q < 7; q++)
					{
						position[0] = cut((int)floor(x_ordinates[j] + 0.5) + p - 3, dimX);
						position[1] = cut((int)floor(y_ordinates[i] + 0.5) + q - 3, dimY);
						distance[0] = position[0] - x_ordinates[j];
						distance[1] = position[1] - y_ordinates[i];

						value2 = bs->compute_gradient_y(distance);
						dif__y += data[position[0] * dimX + position[1]]* value2;

					}
				}
				
				//求二阶导
				for (int p = 0; p < 7; p++)
				{
					for (int q = 0; q < 5; q++)
					{
						position[0] = cut((int)floor(x_ordinates[j] + 0.5) + p - 3, dimX);
						position[1] = cut((int)floor(y_ordinates[i] + 0.5) + q - 3, dimY);
						distance[0] = position[0] - x_ordinates[j];
						distance[1] = position[1] - y_ordinates[i];
						value3 = bs->compute_gradient_xx(distance);
						dif__xx += data[position[0] * dimX + position[1]]* value3;

					}
				}
				for (int p = 0; p < 5; p++)
				{
					for (int q = 0; q < 7; q++)
					{
						position[0] = cut((int)floor(x_ordinates[j] + 0.5) + p - 3, dimX);
						position[1] = cut((int)floor(y_ordinates[i] + 0.5) + q - 3, dimY);
						distance[0] = position[0] - x_ordinates[j];
						distance[1] = position[1] - y_ordinates[i];

						value4 = bs->compute_gradient_yy(distance);
						dif__yy += data[position[0] * dimX + position[1]] * value4;

					}
				}
				for (int p = 0; p < 6; p++)
				{
					for (int q = 0; q < 6; q++)
					{
						position[0] = cut((int)floor(x_ordinates[j] + 0.5) + p - 3, dimX);
						position[1] = cut((int)floor(y_ordinates[i] + 0.5) + q - 3, dimY);
						distance[0] = position[0] - x_ordinates[j];
						distance[1] = position[1] - y_ordinates[i];
						value5 = bs->compute_gradient_xy(distance);

						dif__xy += data[position[0] * dimX + position[1]] * value5;

					}
				}
				/*if (dif__xx > 1)dif__xx = 1;
				if (dif__xx < -1)dif__xx = -1;
				if (dif__yy > 1)dif__yy = 1;
				if (dif__yy < -1)dif__yy = -1;*/
				
				//286_214数据
				/*if(dif__x<3e-14&&dif__x>-3e-14)
				dif__x=0;
				if(dif__y<3e-14&&dif__y>-3e-14)
				dif__y=0;
				if(dif__x==0&&dif__y==0)
				{
				double k1=(dif__xx+dif__yy+sqrt((dif__xx+dif__yy)*(dif__yy+dif__xx)-4*(dif__xx*dif__yy-dif__xy*dif__xy)))/2;
				double k2=(dif__xx+dif__yy-sqrt((dif__xx+dif__yy)*(dif__yy+dif__xx)-4*(dif__xx*dif__yy-dif__xy*dif__xy)))/2;
				if(k1<0&&k1<0)
				Value=0;
				else if(k1>0&&k1>0)
				Value=62.5;
				else Value=115;
				}*/
				//201数据

				if (dif__x< LARGE_ZERO_EPSILON &&dif__x>-LARGE_ZERO_EPSILON)
						dif__x = 0;
				if (dif__y<LARGE_ZERO_EPSILON&&dif__y>-LARGE_ZERO_EPSILON)
						dif__y = 0;
				if (dif__xx<SMALL_ZERO_EPSILON&&dif__xx>-SMALL_ZERO_EPSILON)
					dif__xx = 0;
				if (dif__yy<SMALL_ZERO_EPSILON&&dif__yy>-SMALL_ZERO_EPSILON)
					dif__yy = 0;
				if (dif__xy<SMALL_ZERO_EPSILON&&dif__xy>-SMALL_ZERO_EPSILON)
					dif__xy = 0;

			        CriticalPoint cp;
					cp.xy_local.first = i;
					cp.xy_local.second = j;
					cp.meshIndex = i*dimX + j;
					//cp.dif = make_pair(dif__x*20, dif__y*20);//401sin数据 步长为0.05，偏导数要扩大20倍
					cp.dif = make_pair(dif__x , dif__y );
					
				if (dif__x == 0 && dif__y == 0)
				{
					
					double eig_value[2], eig_vector[2][2], dif[2][2], eps = 1e-50;
					const int Dim = 2, nJt = 1e20;
					
					dif[0][0] = dif__xx; dif[0][1] = dif__xy; dif[1][0] = dif__xy; dif[1][1] = dif__yy;
					//cout << i << " " << j << endl;
					if (JacbiCor(*dif, *eig_vector, eig_value, eps, nJt))
					{
						if (eig_value[1]<0 && eig_value[0]<0)
						{
							cp.type = MAXIMAL;
							cp.eig_vector1 = make_pair(0,0);
							cp.eig_vector2 = make_pair(0,0);
						
							maxPoint.push_back(cp);
						}
						else if (eig_value[0]>0 && eig_value[1]>0)
						{
							cp.type = MINIMAL;
							cp.eig_vector1 = make_pair(0, 0);
							cp.eig_vector2 = make_pair(0, 0);
							
							minPoint.push_back(cp);
						}
						else 
						{
							cp.type = SADDLE;
							//特征向量
							cp.eig_vector1 = make_pair(eig_vector[0][0], eig_vector[1][0]);
							cp.eig_vector2 = make_pair(eig_vector[0][1], eig_vector[1][1]);
							
							saddles.push_back(cp);
						}
						//value <<i << " "<<j<<" "<<cp.type<<endl;
					}	
					keyPoint.push_back(cp);
				}
				else 
					cp.type = REGULAR;

				cp_vec.push_back(cp);
				
				
			
				value << Value << " ";
				/*dif_x << dif__x << " ";
				dif_y << dif__y << " ";
				dif_xx << dif__xx << " ";
				dif_xy << dif__xy << " ";
				dif_yy << dif__yy << " ";*/
			}
			value << endl;
			/*dif_x << endl;
			dif_y << endl;
			dif_xx << endl;
			dif_xy << endl;
			dif_yy << endl;*/
		}
		delete bs;
		cout <<"Critical point number:"<< keyPoint.size() << endl;
		cout << "Max point number:" << maxPoint.size() << endl;
		cout << "Saddle point number:" << saddles.size() << endl;
		cout << "Min point number:" << minPoint.size() << endl;
		cout << "End !" << endl << endl;
	}
	
	bool CPFinder::JacbiCor(double * pMatrix, double *pdblVects, double *pdbEigenValues, double dbEps, int nJt)
	{
		if (pMatrix[0] == pMatrix[3])
		{
			double eig1[2] = { 1, 0 }, eig2[2] = {0, 1};

			if (pMatrix[1]==0&&pMatrix[2]==0)
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

	/*double EvaluateVal(const double x, const double y)
	{
		box_spline* bs;
		bs = new box_spline();
		int position[2];
		double tmpValue, val = 0;
		double distance[2];
		
		if (x<3 || y<3 || x>dimX - 3 || y>dimY - 3)
		{
			return vr.getData(floor(x + 0.5), floor(position[1] + 0.5));
		}
		for (int p = 0; p < 7; p++)
		{
			for (int q = 0; q < 7; q++)
			{
				position[0] = cut((int)floor(x + 0.5) + p - 3, dimX);
				position[1] = cut((int)floor(y + 0.5) + q - 3, dimY);
				distance[0] = position[0] - x;
				distance[1] = position[1] - y;
				tmpValue = bs->compute_value(distance);
				val += vr.getData(position[0], position[1])*tmpValue;
			}
		}
		delete bs;
		if (val > 255)val = 255;
		if (val < 0)val = 0;

		return val;

	}
	double EvaluateGradX(const double x, const double y)
	{
		box_spline* bs;
		bs = new box_spline();
		int position[2];
		double tmpValue, val = 0;;
		double distance[2];

		if ( x<3 ||  y<3 ||  x>dimX - 3 ||  y>dimY - 3)
		{
			return vr.getData(floor( x + 0.5), floor(position[1] + 0.5));
		}
		for (int p = 0; p < 6; p++)
		{
			for (int q = 0; q < 7; q++)
			{
				position[0] = cut((int)floor( x + 0.5) + p - 3, dimX);
				position[1] = cut((int)floor( y + 0.5) + q - 3, dimY);
				distance[0] = position[0] -  x;
				distance[1] = position[1] -  y;
				tmpValue = bs->compute_gradient_x(distance);
				val += vr.getData(position[0], position[1])*tmpValue;
			}
		}
		delete bs;
		if (val > 1)val = 1;
		if (val < -1)val = -1;
		if (val< LARGE_ZERO_EPSILON &&val>-LARGE_ZERO_EPSILON)
			val = 0;
		
		return val;
	}
	double EvaluateGradY(const double x, const double y)
	{
		box_spline* bs;
		bs = new box_spline();
		int position[2];
		double tmpValue, val = 0;;
		double distance[2];

		if ( x<3 ||  y<3 ||  x>dimX - 3 ||  y>dimY - 3)
		{
			return vr.getData(floor( x + 0.5), floor(position[1] + 0.5));
		}
		for (int p = 0; p < 6; p++)
		{
			for (int q = 0; q < 7; q++)
			{
				position[0] = cut((int)floor( x + 0.5) + p - 3, dimX);
				position[1] = cut((int)floor( y + 0.5) + q - 3, dimY);
				distance[0] = position[0] -  x;
				distance[1] = position[1] -  y;
				tmpValue = bs->compute_gradient_y(distance);
				val += vr.getData(position[0], position[1])*tmpValue;
			}
		}
		delete bs;
		if (val > 1)val = 1;
		if (val < -1)val = -1;
		if (val<LARGE_ZERO_EPSILON&&val>-LARGE_ZERO_EPSILON)
			val = 0;
		return val;
	}*/

}