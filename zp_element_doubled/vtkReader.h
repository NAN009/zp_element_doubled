#ifndef _VTK_READER_H
#define _VTK_READER_H
#include "minmax.h"
#include <fstream>
#pragma warning (disable : 4996)

class vtkReader {
public:
	vtkReader() : dataArr(NULL) {}
	~vtkReader() {
		if ( dataArr != NULL ) {
			delete []dataArr;
		}
		if ( upperArr != NULL ) {
			delete []upperArr;
		}
		if ( lowerArr != NULL ) {
			delete []lowerArr;
		}
	}
	const int getSize() const {
		return size;
	}
	const int* getDimension() const {
		return dim;
	}
	const double* getDataArray() const {
		return dataArr;
	}
	const double getData(int x, int y) const {
		return dataArr[ getIndex(x,y) ];
	}
	const double getUpper(int x, int y) const {
		return upperArr[ getIndex(x,y) ];
	}
	const double getLower(int x, int y) const {
		return lowerArr[ getIndex(x,y) ];
	}
	int getIndex(int x, int y) const {
		//先读列
		return max(min(x,dim[0]-1),0) + max(min(y,dim[1]-1),0) * dim[0];
	}//返回（i,j）为数组中第几个数
	bool loadFile(const char* filename);
	int dim[2];
private:
	
	double *dataArr, *upperArr, *lowerArr;
	int size;

	
	
};

bool vtkReader::loadFile(const char* filename) {
	char cachedFilename[1024];
	strcpy(cachedFilename, filename);
	strcat(cachedFilename, ".cache");

	FILE *fp = fopen(cachedFilename, "rb");
	if ( fp!=NULL ) 
	{
		printf("Reading from cache file.\n");
		fread(dim, sizeof(dim), 2, fp);
		printf("Dimension: %d * %d\n", dim[0], dim[1]);
		size = dim[0]*dim[1];
		dataArr = new double[size];
		upperArr = new double[ size ];
		lowerArr = new double[ size ];
		fread(dataArr, sizeof(double), size, fp);
		fread(upperArr, sizeof(double), size, fp);
		fread(lowerArr, sizeof(double), size, fp);
		fclose(fp);
		return true;
	}
	else 
	{
		fp = fopen(filename, "r");
		if ( fp==NULL ) {
			fprintf(stderr, "File %s not found.\n", filename);
			return false;
		}
	}
	char ignored[1024];
	for(int i=0; i<4; i++) {
		fgets(ignored, sizeof(ignored), fp);
		printf("Ignored line: %s", ignored);
	}
	fscanf(fp, "%s%d%d", ignored, &dim[0], &dim[1]);
	//fin >> ignored >> dim[0] >> dim[1] >> dim[2];
	printf("Dimension information: %d * %d\n", dim[0], dim[1]);
	for(int i=0; i<6; i++) {
		fgets(ignored, sizeof(ignored), fp);
		printf("Ignored line: %s", ignored);
	}
	size = dim[0]*dim[1];
	dataArr = new double[ size ];
	upperArr = new double[ size ];
	lowerArr = new double[ size ];
	for(int i=0; i<size; i++) {
		fscanf(fp, "%lf", dataArr+i);
	}
	fclose(fp);
	// Calculating upper and lower bound
	for(int i=0; i<dim[0]; i++) 
	{
		for(int j=0; j<dim[1]; j++) 
		{
			//for(int k=0; k<dim[2]; k++) {
				int index = getIndex(i,j);
				upperArr[index] = -1e300;
				lowerArr[index] = 1e300;
				for(int a=-2; a<=2; a++)
					for(int b=-2; b<=2; b++)
					{
						double data = getData(i+a,j+b);
						if ( upperArr[index] < data )
						{
							upperArr[index] = data;
						}
						if ( lowerArr[index] > data )
						{
							lowerArr[index] = data;
						}
				   }
		}
	}
	// Writing to cache file
	fp = fopen(cachedFilename, "wb");
	fwrite(dim, sizeof(dim), 2, fp);
	fwrite(dataArr, sizeof(double), size, fp);
	fwrite(upperArr, sizeof(double), size, fp);
	fwrite(lowerArr, sizeof(double), size, fp);
	fclose(fp);
	return true;
}

#endif