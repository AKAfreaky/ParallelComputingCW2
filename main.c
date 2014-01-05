#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include "arrayHelpers.h"
#include "mpi.h"


// Global 'debugging' variable
int __VERBOSE;
// John Simms' task ID
#define MASTER 0

// Globals for MPI info
int taskID, numTasks;

#define TAG_INIT_DATA		0
#define TAG_COMPLETE_DATA	1
#define TAG_DATA_CHANGE		2



/* Averages the values surrounding cardinal values in inArray and sets the
	average to outArray. Ignores edges.
	Returns 1 if result changed by less than the precission
*/
int averageFour( double** inArray, double** outArray, int arrayX, int arrayY, double precision)
{
	int i, j, retVal = 1;
	//from pos 1 to arraySize-2 as edges are fixed
	for(i = 1; i < arrayX - 1; i++)
	{
		for(j = 1; j < arrayY - 1; j++)
		{
			double n,s,e,w, oldVal, newVal;

			n	=	inArray[i-1][j];
			s	=	inArray[i+1][j];
			e	=	inArray[i][j+1];
			w	=	inArray[i][j-1];

			outArray[i][j] = (n + s + e + w) / 4.0;

			oldVal = inArray[i][j];
			newVal = outArray[i][j];

			if( fabs(oldVal - newVal) > precision )
			{
				retVal = 0;	// the arrays are too different
			}
		}
	}

	return retVal;
}


int relaxation(double** inArray, double** outArray, int arrayX, int arrayY, double precision)
{
	MPI_Barrier(MPI_COMM_WORLD);// Wait until everyone is good to go.

	int sigStop		= 0;
	int noChange	= 0;

	while( sigStop == 0 )
	{

		noChange = averageFour(inArray, outArray, arrayX, arrayY, precision);

		MPI_Allreduce(&noChange, &sigStop, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

		if (sigStop != 0)
		{
			break;
		}

		if (taskID != MASTER)
			MPI_Send(outArray[1], arrayY, MPI_DOUBLE, taskID - 1, TAG_DATA_CHANGE, MPI_COMM_WORLD);

		if (taskID != (numTasks - 1))
			MPI_Send(outArray[arrayX-2], arrayY, MPI_DOUBLE, taskID + 1, TAG_DATA_CHANGE, MPI_COMM_WORLD);

		if (taskID != MASTER)
			MPI_Recv((void*)outArray[0], arrayY, MPI_DOUBLE, taskID - 1, TAG_DATA_CHANGE, MPI_COMM_WORLD, 0);

		if (taskID != (numTasks - 1))
			MPI_Recv((void*)outArray[arrayX-1], arrayY, MPI_DOUBLE, taskID + 1, TAG_DATA_CHANGE, MPI_COMM_WORLD, 0);

		int i;
		for( i = 0; i < arrayX - 1; i++)
		{
			memcpy(inArray[i], outArray[i], arrayY * sizeof(double));
		}

	}

	return 0;
}


/* Hopefully useful information on how to run the program.
*/
void printUsage()
{
	printf("Arguments are:\n"
			"\t-s\t:\tInteger\t-\tThe size of the matrix\n"
			"\t-p\t:\tdouble\t-\tThe precision to work to\n"
			"\t-r\t:\tInteger\t-\tSeed to use when filling the array. "
			"Zero will use current time() as the seed\n"
			"\t-v\t:\tNone\t-\tFlag to enable more console spew\n"
			"\t-c\t:\tNone\t-\tFlag to enable the correctness test\n");
	// Windows command to stop console applications closing immediately.
	system("pause");
	exit(0);
}


int main(int argc, char **argv)
{
	// Initial values (should get from cmd line)
	int arraySize		= 10;
	double precision	= 10;
	int arrSeed			= time(0);
	int testRight		= 0;

	__VERBOSE 		= 0;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskID);
	printf ("MPI task %d/%d has started...\n", taskID, numTasks);

	// Read options
	// -s is the size, -p is the precision and -t is number of threads.
	// -v turns on some debug spew
	int c;
	opterr = 0;
	while ((c = getopt (argc, argv, "s:p:vr:c")) != -1)
	{
		switch (c)
		{
			case 'p':
            	if (sscanf(optarg, "%lf", &precision) != 1)
				{
					fprintf (stderr,
					"Option -%c requires a double argument.\n",
					 optopt);
					printUsage();
				}
            	break;
			case 'r':
				if (sscanf(optarg, "%i", &arrSeed) != 1)
				{
					fprintf (stderr,
					"Option -%c requires an integer argument.\n",
					 optopt);
					printUsage();
				}
				break;
			case 's':
				if (sscanf(optarg, "%i", &arraySize) != 1)
				{
					fprintf (stderr,
					"Option -%c requires an integer argument.\n",
					 optopt);
					printUsage();
				}
            	break;
			case 'v':
				__VERBOSE = 1;
				break;
		case 'c':
			testRight = 1;
			break;
          	default:
          		printUsage();
           }
	}



	printf("Starting to relax %d square array to precision %f.",
			arraySize, precision);
	printf(" Seed: %d.\n", arrSeed);

	int chunkSize = (arraySize / numTasks) + 2;

	if( taskID == (numTasks - 1) ) //last task
	{
		chunkSize = arraySize - ((arraySize / numTasks) * (numTasks - 1));
	}

	double** currArray;
	double** nextArray;

	// == Setup the data ==
	if (taskID == MASTER)
	{
		currArray = make2DDoubleArray(arraySize, arraySize);
		nextArray = make2DDoubleArray(arraySize, arraySize);
		initArray(currArray, arraySize, arrSeed);

		int	columnPos = chunkSize,  // Skip the chunk the master gets
			i;
		for ( i = 1; i < numTasks; i++) // don't send to self
		{
			if( i == (numTasks - 1) ) //last task
			{
				chunkSize = arraySize - ((arraySize / numTasks) * (numTasks - 1));
			}

			double* oData = malloc(arraySize*chunkSize*sizeof(double));

			int j, currPos = 0;
			for( j = 0; j < chunkSize; j++)
			{
				memcpy(&oData[currPos], currArray[columnPos++], arraySize*sizeof(double));
				currPos += arraySize;
			}

			MPI_Send((void*)oData, (arraySize*chunkSize), MPI_DOUBLE, i, TAG_INIT_DATA, MPI_COMM_WORLD);

			columnPos -= 2; // Overlap the data

			free(oData);
		}

		chunkSize = (arraySize / numTasks) + 2; // we changed the chunkSize val

	}
	else
	{
		double* inBuff = malloc(arraySize*chunkSize*sizeof(double));

		MPI_Recv((void*)inBuff, arraySize*chunkSize, MPI_DOUBLE, MASTER, TAG_INIT_DATA, MPI_COMM_WORLD, 0);

		currArray = make2DDoubleArray(chunkSize, arraySize);
		nextArray = make2DDoubleArray(chunkSize, arraySize);

		int currPos = 0, j;
		for( j = 0; j < chunkSize; j++)
		{
			memcpy(currArray[j], &inBuff[currPos], arraySize*sizeof(double));
			currPos += arraySize;
		}

		free(inBuff);
	}

	// Do the work
	relaxation(currArray, nextArray, chunkSize, arraySize, precision);


	if (taskID != MASTER)
	{
		double* outData = malloc(chunkSize*arraySize*sizeof(double));
		int j, currPos = 0;
		for( j = 2; j < chunkSize; j++) // skip first 2 columns
		{
			memcpy(&outData[currPos], nextArray[j++], arraySize*sizeof(double));
			currPos += arraySize;
		}

		MPI_Send(outData, arraySize*chunkSize, MPI_DOUBLE, MASTER, TAG_COMPLETE_DATA, MPI_COMM_WORLD);
		free(outData);
	}
	else
	{
		int i, columnPos = chunkSize;
		double* inData = malloc(chunkSize*arraySize*sizeof(double));
		for( i = 1; i < numTasks; i++)
		{
			if (i == (numTasks - 1))
			{
				chunkSize = arraySize - ((arraySize / numTasks) * (numTasks - 1));
			}

			MPI_Recv((void*)inData, arraySize*chunkSize, MPI_DOUBLE, i, TAG_COMPLETE_DATA, MPI_COMM_WORLD, 0);

			int currPos = 0, j;
			for( j = 0; j < chunkSize - 2; j++)
			{
				memcpy(nextArray[columnPos++], &inData[currPos], arraySize*sizeof(double));
				currPos += arraySize;
			}
		}
		free(inData);
	}



	// Cleanup
	free2DDoubleArray(currArray, arraySize);
	free2DDoubleArray(nextArray, arraySize);

	double eTime = MPI_Wtime();

	printf("Relaxed %d square matrix in %f seconds\n",
			arraySize, eTime);


	// Windows command to stop console applications closing immediately.
	system("pause");

	return 0;
}
