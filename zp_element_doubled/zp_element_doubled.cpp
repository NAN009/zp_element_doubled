#include <iostream>
#include<fstream>
#include <stdlib.h>

#include "mscomplex.h"

using namespace std;
//void reconstructor::set_dimension(int L,int W)
//{
//	this->L = L;
//	this->W = W;
//}
//
//void reconstructor::set_step(double stepx,double stepy)
//{
//	this->stepx = stepx;
//	this->stepy = stepy;
//
//    this->nx = (int)floor(((this->L)-1)/this->stepx)+1;
//	this->ny = (int)floor(((this->W)-1)/this->stepy)+1;
//}
//
//int cut(int i,int N)
//{
//	int tmp = ((i<0)?0:i);
//	return ((tmp>(N-1))?(N-1):tmp);
//}


int main()
{
	msc2d::MSComplex2D msc;
	msc.createMSComplex2D();
	msc.saveMSComplex();
	//box_spline* bs;
	//bs = new box_spline();
	////FILE * fp;
	//
	//vr.loadFile("D:\\data\\201sin.vtk");//201数据：由于201数据步长为0.1，求偏导数的默认步长为1，故求出结果为实际结果的十分之一；
	//double x_ordinates[10000],y_ordinates[10000];//第一个参数为列，按列读取数据
	//ofstream value("D:\\data\\201sin_value.txt");
	//ofstream dif_x("D:\\data\\201sin_dif_x.txt");
	//ofstream dif_y("D:\\data\\201sin_dif_y.txt");
	//ofstream dif_xx("D:\\data\\201sin_dif_xx.txt");
	//ofstream dif_xy("D:\\data\\201sin_dif_xy.txt");
	//ofstream dif_yy("D:\\data\\201sin_dif_yy.txt");
	////cout<<vr.getData(0,0)<<" "<<vr.getData(0,1)<<" "<<vr.getData(0,2)<<endl;
	////cout<<vr.getData(1,0)<<" "<<vr.getData(1,1)<<" "<<vr.getData(1,2)<<endl;
	//int m=0;
	//int position[2];
	//double sum,sum1,sum2,sum4,sum3,sum5,sum0;
	//double distance[2];
	//double diff_x,diff_xx;
	//double diff_y,diff_yy;
	//double f_with_diff,f_with_diff2;
	//double value1,value3,value2,value4,value5,value0;
	//for(int i=0;i<vr.dim[0];i++)
	//	x_ordinates[i] = i;
	//for(int i=0;i<vr.dim[1];i++)
	//	y_ordinates[i] = i;
	//
	//for(int i=0;i<vr.dim[0];i++)	
	//{
	//	for(int j=0;j<vr.dim[1];j++)
	//	{
	//		sum=0;sum1=0;sum2=0;sum4=0,sum5=0,sum0=0,sum3=0;
	//		//value
	//		if(i==0||i==1||j==0||j==1||i==vr.dim[0]-1||i==vr.dim[0]-2||j==vr.dim[1]-1||j==vr.dim[1]-2)
	//			sum0 = vr.getData(position[0],position[1]);
	//		else
	//		{
	//			for(int p=0;p<7;p++)
	//			{
	//				for(int q=0;q<7;q++)
	//				{
	//					position[0]=cut((int)floor(x_ordinates[i]+0.5)+p-3,vr.dim[0]);
	//					position[1] = cut((int)floor(y_ordinates[j]+0.5)+q-3,vr.dim[1]);
	//					distance[0] = position[0]-x_ordinates[i];
	//					distance[1] = position[1]-y_ordinates[j];
	//					value0 = bs->compute_value(distance);
	//				//cout<<p7.evaluate(distance)<<endl;
	//				//diff_x = vr.getData(cut(position[0]+1,vr.dim[0]),position[1])+vr.getData(cut(position[0]-1,vr.dim[0]),position[1])-2*vr.getData(position[0],position[1]);
	//				//diff_y = vr.getData(position[0],cut(position[1]+1,vr.dim[1]))+vr.getData(position[0],cut(position[1]-1,vr.dim[1]))-2*vr.getData(position[0],position[1]);
	//				//f_with_diff = vr.getData(position[0],position[1])-5.0*(diff_x+diff_y)/24;
	//				//diff_x = vr.getData(cut(position[0]+1,vr.dim[0]),position[1])+vr.getData(cut(position[0]-1,vr.dim[0]),position[1])-2*vr.getData(position[0],position[1]);
	//				//diff_y = vr.getData(position[0],cut(position[1]+1,vr.dim[1]))+vr.getData(position[0],cut(position[1]-1,vr.dim[1]))-2*vr.getData(position[0],position[1]);
	//				//f_with_diff = vr.getData(position[0],position[1])-1.0*(diff_x+diff_y)/8;
	//				
	//					sum0 += vr.getData(position[0],position[1])*value0;
	//				}
	//			}
	//		}
	//		if (sum0>4)sum0=4;
	//		//求导
	//		for(int p=0;p<6;p++)
	//		{
	//			for(int q=0;q<7;q++)
	//			{
	//				position[0]=cut((int)floor(x_ordinates[i]+0.5)+p-3,vr.dim[0]);
	//				position[1] = cut((int)floor(y_ordinates[j]+0.5)+q-3,vr.dim[1]);
	//				distance[0] = position[0]-x_ordinates[i];
	//				distance[1] = position[1]-y_ordinates[j];
	//					
	//				value1 = bs->compute_gradient_x(distance);
	//				value2 = bs->compute_gradient_y(distance);
	//				diff_x = vr.getData(cut(position[0]+1,vr.dim[0]),position[1])+vr.getData(cut(position[0]-1,vr.dim[0]),position[1])-2*vr.getData(position[0],position[1]);
	//				diff_y = vr.getData(position[0],cut(position[1]+1,vr.dim[1]))+vr.getData(position[0],cut(position[1]-1,vr.dim[1]))-2*vr.getData(position[0],position[1]);
	//				f_with_diff = vr.getData(position[0],position[1])-1.0*(diff_x+diff_y)/8;
	//				sum1 += f_with_diff*value1;
	//				sum2 += f_with_diff*value2;
	//				//sum3 += vr.getData(position[0],position[1])*value3;
	//				//sum4 += vr.getData(position[0],position[1])*value4;
	//			}
	//		}
	//		//导数最值
	//		if (sum1 > 1)sum1 = 1;
	//		if (sum1 < -1)sum1 = -1;
	//		if (sum2 > 1)sum2 = 1;
	//		if (sum2< -1)sum2 = -1;
	//		//求二阶导
	//		for(int p=0;p<7;p++)
	//		{
	//			for(int q=0;q<5;q++)
	//			{
	//				position[0]=cut((int)floor(x_ordinates[i]+0.5)+p-3,vr.dim[0]);
	//				position[1] = cut((int)floor(y_ordinates[j]+0.5)+q-3,vr.dim[1]);
	//				distance[0] = position[0]-x_ordinates[i];
	//				distance[1] = position[1]-y_ordinates[j];
	//				value3 = bs->compute_gradient_xx(distance);	
	//			//	value4 = bs->compute_gradient_yy(distance);
	//				//diff_xx = vr.getData(cut(position[0]+1,vr.dim[0]),position[1])+vr.getData(cut(position[0]-1,vr.dim[0]),position[1])-2*vr.getData(position[0],position[1]);
	//				//diff_yy = vr.getData(position[0],cut(position[1]+1,vr.dim[1]))+vr.getData(position[0],cut(position[1]-1,vr.dim[1]))-2*vr.getData(position[0],position[1]);
	//				//f_with_diff2 = vr.getData(position[0],position[1])-1.0*(diff_x+diff_y)/8;
	//				//sum3 += f_with_diff2*value3;
	//				//sum4 += f_with_diff2*value4;
	//				sum3 += vr.getData(position[0],position[1])*value3;
	//				//sum4 += vr.getData(position[0],position[1])*value4;
	//				
	//			}
	//		}
	//		for(int p=0;p<5;p++)
	//		{
	//			for(int q=0;q<7;q++)
	//			{
	//				position[0]=cut((int)floor(x_ordinates[i]+0.5)+p-3,vr.dim[0]);
	//				position[1] = cut((int)floor(y_ordinates[j]+0.5)+q-3,vr.dim[1]);
	//				distance[0] = position[0]-x_ordinates[i];
	//				distance[1] = position[1]-y_ordinates[j];
	//				//value3 = bs->compute_gradient_xx(distance);	
	//				value4 = bs->compute_gradient_yy(distance);
	//				//diff_xx = vr.getData(cut(position[0]+1,vr.dim[0]),position[1])+vr.getData(cut(position[0]-1,vr.dim[0]),position[1])-2*vr.getData(position[0],position[1]);
	//				//diff_yy = vr.getData(position[0],cut(position[1]+1,vr.dim[1]))+vr.getData(position[0],cut(position[1]-1,vr.dim[1]))-2*vr.getData(position[0],position[1]);
	//				//f_with_diff2 = vr.getData(position[0],position[1])-1.0*(diff_x+diff_y)/8;
	//				//sum3 += f_with_diff2*value3;
	//				//sum4 += f_with_diff2*value4;
	//				//sum3 += vr.getData(position[0],position[1])*value3;
	//				sum4 += vr.getData(position[0],position[1])*value4;
	//				
	//			}
	//		}
	//		for(int p=0;p<6;p++)
	//		{
	//			for(int q=0;q<6;q++)
	//			{
	//				position[0]=cut((int)floor(x_ordinates[i]+0.5)+p-3,vr.dim[0]);
	//				position[1] = cut((int)floor(y_ordinates[j]+0.5)+q-3,vr.dim[1]);
	//				distance[0] = position[0]-x_ordinates[i];
	//				distance[1] = position[1]-y_ordinates[j];
	//				value5 = bs->compute_gradient_xy(distance);
	//				
	//				sum5 += vr.getData(position[0],position[1])*value5;
	//				
	//			}
	//		}
	//		if (sum3 > 1)sum3 = 1;
	//		if (sum3 < -1)sum3 = -1;
	//		if (sum4 > 1)sum4 = 1;
	//		if (sum4 < -1)sum4 = -1;
	//		sum5 = 0;
	//		//286_214数据
	//		/*if(sum1<3e-14&&sum1>-3e-14)
	//			sum1=0;
	//		if(sum2<3e-14&&sum2>-3e-14)
	//			sum2=0;
	//		if(sum1==0&&sum2==0)
	//		{
	//		double k1=(sum3+sum4+sqrt((sum3+sum4)*(sum4+sum3)-4*(sum3*sum4-sum5*sum5)))/2;
	//		double k2=(sum3+sum4-sqrt((sum3+sum4)*(sum4+sum3)-4*(sum3*sum4-sum5*sum5)))/2;
	//		if(k1<0&&k1<0)
	//			sum0=0;
	//		else if(k1>0&&k1>0)
	//			sum0=62.5;
	//		else sum0=115;
	//		}*/
	//		//201数据
	//		
	//		if(sum1<0.006&&sum1>-0.006)
	//			sum1=0;
	//		if(sum2<0.006&&sum2>-0.006)
	//			sum2=0;
	//		if(sum3<2e-006&&sum3>-2e-006)
	//			sum3=0;
	//		if(sum4<2e-006&&sum4>-2e-006)
	//			sum4=0;
	//		if(sum5<2e-006&&sum5>-2e-006)
	//			sum5=0;
	//		//double eigvalue1,eigvalue2;
	//		//if(sum1==0&&sum2==0)
	//		//{
	//		//	eigvalue1=0,eigvalue2=0;
	//		//	
	//		//	eigvalue1=(sum3+sum4+sqrt((sum3+sum4)*(sum4+sum3)-4*(sum3*sum4-sum5*sum5)))/2;
	//		//	eigvalue2=(sum3+sum4-sqrt((sum3+sum4)*(sum4+sum3)-4*(sum3*sum4-sum5*sum5)))/2;
	//		//	
	//		//	if (eigvalue1<1e-10&&eigvalue1>-1e-10)eigvalue1 = 0;
	//		//	if (eigvalue2<1e-10&&eigvalue2>-1e-10)eigvalue2 = 0;
	//		//	if (eigvalue1<0 && eigvalue2<0)
	//		//		sum0=2;		//极大值点
	//		//	else if (eigvalue1>0 && eigvalue2>0)
	//		//		sum0=4;		//极小值点
	//		//	else sum0=0;		//鞍点
	//		//}
	//		
	//		value << sum0<<" ";
	//		dif_x << sum1 << " ";
	//		dif_y << sum2 << " ";
	//		dif_xx << sum3 << " ";
	//		dif_xy << sum5 << " ";
	//		dif_yy << sum4 << " ";
	//		//os2<<sum2<<" ";
	//		//os3<<sum5<<" ";
	//		
	//	}
	//	value << endl;
	//	dif_x << endl;
	//	dif_y << endl;
	//	dif_xx << endl;
	//	dif_xy << endl;
	//	dif_yy << endl;
	//}
	////x,y取值范围都为【-3,3】
	//*box_spline* bs;
	//bs = new box_spline();
	//double a[2];
	//a[0] = 2;
	//a[1] = 0;
	//double value = bs->compute_value(a);
	//double value1 = bs->compute_gradient_x(a);
	//double value2 = bs->compute_gradient_y(a);
	//double value3 = bs->compute_gradient_xx(a);
	//double value4 = bs->compute_gradient_xy(a);
	//double value5 = bs->compute_gradient_yy(a);
	//cout<<value<<endl;
	//cout<<value1<<endl;
	//cout<<value2<<endl;
	//cout<<value3<<endl;
	//cout<<value4<<endl;
	//cout<<value5<<endl;
	//delete bs;*/
	//delete bs;
	system("pause");
	return 0;
}
