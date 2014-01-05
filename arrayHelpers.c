#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// fills a 2D-square array with random values
void initArray( double** theArray, int arraySize, int seed )
{
	srand( seed ? seed : time( 0 ) );
	int i, j;
	for( i = 0; i < arraySize; i++ )
	{
		for( j = 0; j < arraySize; j++ )
		{
			double value =  (rand() % 10000) + 1.0;
			theArray[i][j] = value;
		}
	}
}

void printSquareArray( double** theArray, int arraySize )
{
	int i, j;
	for( i = 0; i < arraySize; i++ )
	{
		for( j = 0; j < arraySize; j++ )
		{
			printf("%f\t", theArray[i][j]);
		}
		printf("\n");
	}
}


double** make2DDoubleArray(int arraySizeX, int arraySizeY)
{
	double **theArray;
	int i;
	theArray = malloc(arraySizeX*sizeof(double*));
	for (i = 0; i < arraySizeX; i++)
	{
		theArray[i] = malloc(arraySizeY*sizeof(double));
	}
	return theArray;
}

void free2DDoubleArray(double** theArray, int arraySizeX)
{
	int i;
	for (i = 0; i < arraySizeX; i++)
	{
   		free(theArray[i]);
	}
	free(theArray);
}
