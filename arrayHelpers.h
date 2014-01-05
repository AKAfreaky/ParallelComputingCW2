#ifndef ARRAYHELPERS_H_INCLUDED
#define ARRAYHELPERS_H_INCLUDED

// fills a 2D-square array with random values
void initArray( double** theArray, int arraySize, int seed );

void initArrayPattern( double** theArray, int arraySize );

void printSquareArray( double** theArray, int arraySize );

double** make2DDoubleArray(int arraySizeX, int arraySizeY);

void free2DDoubleArray(double** theArray, int arraySizeX);

#endif // ARRAYHELPERS_H_INCLUDED
