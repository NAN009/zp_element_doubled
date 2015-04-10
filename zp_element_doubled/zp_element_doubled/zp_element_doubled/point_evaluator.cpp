#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include "point_evaluator.h"

using namespace std;

point_evaluator::point_evaluator(void)
{
}

point_evaluator::~point_evaluator(void)
{
	delete gf;
}

void point_evaluator::initialize(void)
{
	gf = new global_functions();

	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<4;k++)
				inv_V[i][j][k] = 0;

	inv_V[0][0][0] = -1;
	inv_V[0][1][0] = 1;
	inv_V[1][0][0] = -1;
	inv_V[1][1][0] = -1;
	inv_V[1][2][0] = 1;
	inv_V[2][0][0] = 2;
	inv_V[0][0][1] = -1;
	inv_V[0][1][1] = -1;
	inv_V[0][2][1] = 1;
	inv_V[1][0][1] = 1;
	inv_V[1][1][1] = -1;
	inv_V[2][1][1] = 2;
	inv_V[0][0][2] = 1;
	inv_V[0][1][2] = 1;
	inv_V[0][2][2] = -1;
	inv_V[1][0][2] = -1;
	inv_V[1][1][2] = 1;
	inv_V[2][1][2] = -2;
	inv_V[2][2][2] = 2;
	inv_V[0][0][3] = 1;
	inv_V[0][1][3] = -1;
	inv_V[1][0][3] = 1;
	inv_V[1][1][3] = 1;
	inv_V[1][2][3] = -1;
	inv_V[2][0][3] = -2;
	inv_V[2][2][3] = 2;

    ifstream is;
	char* file_name = "ctrl_points.txt";
	try
	{
	    is.open(file_name,fstream::in);//以输入的方式打开文件
	}
	catch(exception e)
	{
		std::cerr<<"error opening the file"<<std::endl;
		system("pause");
		exit(1);
	}

	string line;
	char* cell;
	char* num;
	char* den;

	char* outer = NULL;
	char* inner = NULL;

	for(int i=0;i<6;i++)
	{
		for(int j=0;j<6;j++)
		{
			for(int row_num=0;row_num<28;row_num++)
			{	
				getline(is,line);//将is 的内容读取到line 中，遇到换行符结束
				outer = NULL;
				cell = strtok_s(&(line[0])," ",&outer);//分解字符串为一组字符串；分解line[0],并用空格隔开，将剩余的存入outer,得到cell为cp的具体值
				for(int col_num=0;col_num<4;col_num++)
				{
					if(NULL != cell)
					{
						if(0 != strcmp(cell,"0"))//if cell值为不为0时，执行
						{
							num = strtok_s(cell,"/",&inner); //以下两行实现 分割cell
							den = strtok_s(NULL,"/",&inner);
							this->cp[0][row_num][col_num][i][j] = gf->str2int((string)num);//将cell 的值赋给cp[0],cp[1]存放“/0”
							this->cp[1][row_num][col_num][i][j] = gf->str2int((string)den);
						}
						else
						{
							this->cp[0][row_num][col_num][i][j] = 0;
							this->cp[1][row_num][col_num][i][j] = 1;
						}
					}
					else
					{
						cerr<<"error,column number less than 24!"<<endl;
						system("pause");
						exit(1);
					}
					cell = strtok_s(NULL," ",&outer);
				}
			}
		}
	}

	is.close();

	Ne[0][0] = 1;
	Ne[0][1] = -1;
	Ne[1][0] = 1;
	Ne[1][1] = 1;

}

double point_evaluator::evaluate(double *y)
{
  double *tmp;
  tmp = y;
  if(out_boundary(tmp))
    return 0.0;

  this->x[0] = *y;
  y++;
  this->x[1] = *y;

  int i_square[2];
  double x_square[2];

  for(int i=0;i<2;i++)
  {
    i_square[i] = (int)floor(x[i]);//（int）向下取整，int(-.39)=-4
    x_square[i] = x[i] - i_square[i];
  }

  double values[2];
  int index[2];

  values[0] = Ne[0][0]*x_square[0]+Ne[0][1]*x_square[1];//相减
  values[1] = Ne[1][0]*x_square[0]+Ne[1][1]*x_square[1]-1;//相加-1

  index[0] = (values[0]>=0)?1:0;//values为正数或0时 index为1 否则为0 
  index[1] = (values[1]>=0)?1:0;

  int triangle_index;

  if(index[0]==0)
  {
	  if(index[1]==0)
		  triangle_index=0;
	  else
		  triangle_index = 2;
  }
  else
  {
	  if(index[1]==0)
		  triangle_index = 1;
	  else
		  triangle_index = 3;
  }

  //bary为三维矩阵中t_i列的值之和
  for(int i=0;i<3;i++)
  {
    bary[i] = inv_V[i][0][triangle_index]*x_square[0]+inv_V[i][1][triangle_index]*x_square[1]+inv_V[i][2][triangle_index];
  }

  double value = 0;
  int n=0;
  int k;

  double temp1; //temp;
  double temp2;
  double temp3;
  double temp4;
  double temp5;
  double temp;

  for(int i=0;i<=6;i++)
  {
    for(int j=0;j<=(6-i);j++)
    {
      k = 6-i-j;

	  temp1 = cp[0][n][triangle_index][i_square[0]+3][i_square[1]+3];
      temp2 = cp[1][n][triangle_index][i_square[0]+3][i_square[1]+3];
	  temp3 = my_pow(bary[1],j);
	  temp4 = my_pow(bary[0],i);
	  temp5 = my_pow(bary[2],k);
	  temp = temp1*720.0/(factorial(i)*factorial(j)*factorial(k)*temp2)*temp3*temp4*temp5;// 6!=720 double_zp为6次，所以取6 
      value += temp;

      n++;
    }
  }

  return value;
}

int point_evaluator::factorial(int i)
{
  if(i>0)
    return i*factorial(i-1);
  else
    return 1;
}

bool point_evaluator::out_boundary(double* tmp)
{
  bool flag = false;

  if(*tmp<=-3 || *tmp>=3)
    flag = true;
  tmp++;
  if(*tmp<=-3 || *tmp>=3)
    flag = true;

  return flag;
}

double point_evaluator::my_pow(double base,int power)
{
	if(power==0)
		return 1.0;
	else
		return base*my_pow(base,power-1);
}