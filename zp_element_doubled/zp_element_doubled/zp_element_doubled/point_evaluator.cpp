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
	    is.open(file_name,fstream::in);//������ķ�ʽ���ļ�
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
				getline(is,line);//��is �����ݶ�ȡ��line �У��������з�����
				outer = NULL;
				cell = strtok_s(&(line[0])," ",&outer);//�ֽ��ַ���Ϊһ���ַ������ֽ�line[0],���ÿո��������ʣ��Ĵ���outer,�õ�cellΪcp�ľ���ֵ
				for(int col_num=0;col_num<4;col_num++)
				{
					if(NULL != cell)
					{
						if(0 != strcmp(cell,"0"))//if cellֵΪ��Ϊ0ʱ��ִ��
						{
							num = strtok_s(cell,"/",&inner); //��������ʵ�� �ָ�cell
							den = strtok_s(NULL,"/",&inner);
							this->cp[0][row_num][col_num][i][j] = gf->str2int((string)num);//��cell ��ֵ����cp[0],cp[1]��š�/0��
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
    i_square[i] = (int)floor(x[i]);//��int������ȡ����int(-.39)=-4
    x_square[i] = x[i] - i_square[i];
  }

  double values[2];
  int index[2];

  values[0] = Ne[0][0]*x_square[0]+Ne[0][1]*x_square[1];//���
  values[1] = Ne[1][0]*x_square[0]+Ne[1][1]*x_square[1]-1;//���-1

  index[0] = (values[0]>=0)?1:0;//valuesΪ������0ʱ indexΪ1 ����Ϊ0 
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

  //baryΪ��ά������t_i�е�ֵ֮��
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
	  temp = temp1*720.0/(factorial(i)*factorial(j)*factorial(k)*temp2)*temp3*temp4*temp5;// 6!=720 double_zpΪ6�Σ�����ȡ6 
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