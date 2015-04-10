#include <iostream>
#include <stdlib.h>

#include "box_spline.h"
#include "point_evaluator_6direction2.h"

using namespace std;

int main()
{
	box_spline* bs;

	bs = new box_spline();

	double a[2];

	
	a[0] = 2;
	a[1] = 1;

	double value = bs->compute_value(a);
	double value1 = bs->compute_gradient_x(a);
	double value2 = bs->compute_gradient_y(a);
	double value3 = bs->compute_gradient_xx(a);
	double value4 = bs->compute_gradient_xy(a);
	double value5 = bs->compute_gradient_yy(a);

	cout<<value<<endl;
	cout<<value1<<endl;
	cout<<value2<<endl;
	cout<<value3<<endl;
	cout<<value4<<endl;
	cout<<value5<<endl;
	

	delete bs;
	
	system("pause");
	
}
