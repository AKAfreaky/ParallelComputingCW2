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
int taskID, numTasks, normChunkSize, lastChunkSize;

#define TAG_INIT_DATA		0
#define TAG_COMPLETE_DATA	1
#define TAG_DATA_CHANGE		2



int checkDiff( double** oldArray, double** newArray, int arrayX, int arrayY, double precision)
{
	int i,j;

	for(i = 1; i < arrayX - 1; i++)
	{
		for(j = 1; j < arrayY - 1; j++)
		{
			double oldVal = oldArray[i][j],
				   newVal = newArray[i][j];

			// if the values differ by more than the precision
			if( fabs(oldVal - newVal) > precision )
			{
				return 0;	// the arrays are too different
			}
		}
	}

	// Assume we passed
	return 1;
}



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
		if (__VERBOSE)
			printf("Task %d starting averaging. maxX:%d, maxY:%d\n", taskID, arrayX, arrayY);

		noChange = averageFour(inArray, outArray, arrayX, arrayY, precision);

		MPI_Allreduce(&noChange, &sigStop, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

		if (sigStop != 0)
		{
			break;
		}

		MPI_Status status;

		if (taskID != MASTER)
			MPI_Send(outArray[1], arrayY, MPI_DOUBLE, taskID - 1, TAG_DATA_CHANGE, MPI_COMM_WORLD);

		if (taskID != (numTasks - 1))
			MPI_Send(outArray[arrayX-2], arrayY, MPI_DOUBLE, taskID + 1, TAG_DATA_CHANGE, MPI_COMM_WORLD);

		if (taskID != MASTER)
			MPI_Recv(outArray[0], arrayY, MPI_DOUBLE, taskID - 1, TAG_DATA_CHANGE, MPI_COMM_WORLD, &status);

		if (taskID != (numTasks - 1))
			MPI_Recv(outArray[arrayX-1], arrayY, MPI_DOUBLE, taskID + 1, TAG_DATA_CHANGE, MPI_COMM_WORLD, &status);

		int i;
		for( i = 0; i < arrayX; i++)
		{
			memcpy(inArray[i], outArray[i], arrayY * sizeof(double));
		}

		MPI_Barrier(MPI_COMM_WORLD);

	}

	return 0;
}


int getChunkSize(int taskRank)
{
	if( taskRank == (numTasks - 1) )
	{
		return lastChunkSize;
	}

	return normChunkSize;
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
	exit(0);
}


int main(int argc, char **argv)
{
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskID);

	// Initial values (should get from cmd line)
	int arraySize		= 100;
	double precision	= 1.5;
	int testRight		= 0;

	__VERBOSE 		= 0;

	// Read options
	// -s is the size, -p is the precision and -t is number of threads.
	// -v turns on some debug spew
	int c;
	opterr = 0;
	while ((c = getopt (argc, argv, "s:p:vc")) != -1)
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

	normChunkSize	= (arraySize / numTasks) + 2;
	lastChunkSize	= arraySize - ((arraySize / numTasks) * (numTasks - 1));
	int chunkSize 	= getChunkSize(taskID);

	double sTime 	= MPI_Wtime();

	double** startArray;
	double** endArray;

	if(__VERBOSE)
	{
		printf("MPI task %d/%d has started...\n", taskID, numTasks);
	}

	// == Setup the data ==
	if (taskID == MASTER)
	{
		if (__VERBOSE)
		{
			printf("Starting to relax %d square array to precision %f.\n",
														arraySize, precision);
		}

		// malloc arrays
		startArray = make2DDoubleArray(arraySize, arraySize);
		endArray = make2DDoubleArray(arraySize, arraySize);

		// fill arrays
		initArrayPattern(startArray, arraySize);
		initArrayPattern(endArray, arraySize);

		if (__VERBOSE)
		{
			printf("Initial array:\n");
			printSquareArray(startArray, arraySize);
		}

		int	columnPos = chunkSize - 2, i;  // Skip the chunk the master gets

		for ( i = 1; i < numTasks; i++) // don't send to self
		{
			chunkSize = getChunkSize(i);

			double* sendBuff = malloc(arraySize*chunkSize*sizeof(double));

			int j, currPos = 0;
			for( j = 0; j < chunkSize; j++)
			{
				memcpy(&sendBuff[currPos], startArray[columnPos++], arraySize*sizeof(double));
				currPos += arraySize;
			}

			MPI_Send(sendBuff, (arraySize*chunkSize), MPI_DOUBLE, i, TAG_INIT_DATA, MPI_COMM_WORLD);

			if(__VERBOSE)
				printf("Master sent data to task %d\n", i);

			columnPos -= 2; // Overlap the data

			free(sendBuff);
		}

		chunkSize = getChunkSize(taskID);// we changed the chunkSize val

	}
	else
	{
		double* inBuff = malloc(arraySize*chunkSize*sizeof(double));
		MPI_Status status;

		if(__VERBOSE)
			printf("Task %d, waiting for data...\n", taskID);

		MPI_Recv((void*)inBuff, arraySize*chunkSize, MPI_DOUBLE, MASTER, TAG_INIT_DATA, MPI_COMM_WORLD, &status);

		if(__VERBOSE)
			printf("Task %d, recieved data...\n", taskID);

		startArray	= make2DDoubleArray(chunkSize, arraySize);
		endArray	= make2DDoubleArray(chunkSize, arraySize);

		int currPos = 0, j;
		for( j = 0; j < chunkSize; j++)
		{
			memcpy(startArray[j], &inBuff[currPos], arraySize*sizeof(double));
			memcpy(endArray[j], startArray[j],	arraySize*sizeof(double));
			currPos += arraySize;
		}

		free(inBuff);
	}

	if(__VERBOSE)
		printf("Starting relaxation for task %d.\n", taskID);

	// == Do the work ==
	relaxation(startArray, endArray, chunkSize, arraySize, precision);

	if(__VERBOSE)
		printf("Finished relaxation for task %d.\n", taskID);

	if (taskID != MASTER)
	{
		double* outData = malloc(chunkSize*arraySize*sizeof(double));

		int j, currPos = 0;
		for( j = 2; j < chunkSize; j++) // skip first 2 columns
		{
			memcpy(&outData[currPos], endArray[j], arraySize*sizeof(double));
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
			if(__VERBOSE)
				printf("Receiving from task %d\n", i);

			// Resize recv buffer if needed
			int oldChunkSize = chunkSize;
			chunkSize = getChunkSize(i);

			if (chunkSize != oldChunkSize)
			{
				inData = realloc(inData, chunkSize*arraySize*sizeof(double));
			}


			MPI_Status status;
			MPI_Recv(inData, arraySize*chunkSize, MPI_DOUBLE, i, TAG_COMPLETE_DATA, MPI_COMM_WORLD, &status);

			int currPos = 0, j;
			for( j = 0; j < chunkSize - 2; j++)
			{
				memcpy(endArray[columnPos++], &inData[currPos], arraySize*sizeof(double));
				currPos += arraySize;
			}
		}
		free(inData);
	}


	if (__VERBOSE && (taskID == MASTER))
	{
		printf("Parallel result:\n");
		printSquareArray(endArray, arraySize);
	}

	if (testRight && (taskID == MASTER))
	{
		double** checkArray1	= make2DDoubleArray(arraySize, arraySize);
		double** checkArray2	= make2DDoubleArray(arraySize, arraySize);

		initArrayPattern(checkArray1, arraySize);
		initArrayPattern(checkArray2, arraySize);

		int fin = 0;

		while( fin == 0 )
		{
			fin = averageFour(checkArray1, checkArray2, arraySize, arraySize, precision);

			int i;
			for( i = 1; i < arraySize - 1; i++)
			{
				memcpy(checkArray1[i], checkArray2[i], arraySize * sizeof(double));
			}
		}

		int diff = checkDiff(endArray, checkArray2, arraySize, arraySize, precision);

		printf("Parallel result %s the serial result\n", diff ? "matched" : "didn't match");
	}

	if(__VERBOSE)
		printf("Starting cleanup for task %d.\n", taskID);

	// Cleanup
	chunkSize = getChunkSize(taskID);

	free2DDoubleArray(startArray, (taskID == MASTER) ? arraySize : chunkSize);
	free2DDoubleArray(endArray, (taskID == MASTER) ? arraySize : chunkSize);

	double eTime = MPI_Wtime();

	printf("Relaxed %d square matrix in %f seconds (task: %d/%d)\n",
			arraySize, (eTime - sTime), taskID, numTasks);


	MPI_Finalize();

	system("pause");

	return 0;
}
