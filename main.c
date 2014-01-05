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
			MPI_Recv((void*)outArray[0], arrayY, MPI_DOUBLE, taskID - 1, TAG_DATA_CHANGE, MPI_COMM_WORLD, &status);

		if (taskID != (numTasks - 1))
			MPI_Recv((void*)outArray[arrayX-1], arrayY, MPI_DOUBLE, taskID + 1, TAG_DATA_CHANGE, MPI_COMM_WORLD, &status);

		int i;
		for( i = 0; i < arrayX - 1; i++)
		{
			memcpy(inArray[i], outArray[i], arrayY * sizeof(double));
		}

	}

	return 0;
}


int calculateChunkSize(int arraySize, int taskRank)
{
	int chunkSize = (arraySize / numTasks) + 2;

	if( taskRank == (numTasks - 1) ) //last task
	{
		chunkSize = arraySize - ((arraySize / numTasks) * (numTasks - 1));
	}

	return chunkSize;
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
	printf ("MPI task %d/%d has started...\n", taskID, numTasks);

	// Initial values (should get from cmd line)
	int arraySize		= 10;
	double precision	= 10;
	int arrSeed			= time(0);
	int testRight		= 0;

	__VERBOSE 		= 0;

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

	int chunkSize = calculateChunkSize(arraySize, taskID);

	double** currArray;
	double** nextArray;

	// == Setup the data ==
	if (taskID == MASTER)
	{
		printf("Starting to relax %d square array to precision %f.\n",
			arraySize, precision);

		printf("Master creating its arrays.\n");
		currArray = make2DDoubleArray(arraySize, arraySize);
		nextArray = make2DDoubleArray(arraySize, arraySize);
		initArrayPattern(currArray, arraySize);
		printf("Master created its arrays.\n");

		if (__VERBOSE)
		{
			printf("Initial array:\n");
			printSquareArray(currArray, arraySize);
		}

		int	columnPos = chunkSize - 2,  // Skip the chunk the master gets
			i;
		for ( i = 1; i < numTasks; i++) // don't send to self
		{
			chunkSize = calculateChunkSize(arraySize, i);

			double* oData = malloc(arraySize*chunkSize*sizeof(double));

			if (oData == NULL)
			{
				printf("malloc error! crashing...");
			}

			int j, currPos = 0;
			for( j = 0; j < chunkSize; j++)
			{
				memcpy(&oData[currPos], currArray[columnPos++], arraySize*sizeof(double));
				currPos += arraySize;
				if (currPos >= arraySize*chunkSize)
				{
					printf("OOps we've copied too much!\n");

				}
			}

			printf("Sending data to task %d\n", i);
			MPI_Send(oData, (arraySize*chunkSize), MPI_DOUBLE, i, TAG_INIT_DATA, MPI_COMM_WORLD);
			printf("Sent data to task %d\n", i);

			columnPos -= 2; // Overlap the data

			free(oData);
		}

		chunkSize = calculateChunkSize(arraySize, taskID);// we changed the chunkSize val

	}
	else
	{
		double* inBuff = malloc(arraySize*chunkSize*sizeof(double));
		MPI_Status status;

		printf("Task %d, waiting for data...\n", taskID);
		MPI_Recv((void*)inBuff, arraySize*chunkSize, MPI_DOUBLE, MASTER, TAG_INIT_DATA, MPI_COMM_WORLD, &status);
		printf("Task %d, recieved data...\n", taskID);

		currArray = make2DDoubleArray(chunkSize, arraySize);
		nextArray = make2DDoubleArray(chunkSize, arraySize);

		int currPos = 0, j;
		for( j = 0; j < chunkSize; j++)
		{
			memcpy(currArray[j], &inBuff[currPos], arraySize*sizeof(double));
			memcpy(nextArray[j], currArray[j],	arraySize*sizeof(double));
			currPos += arraySize;
		}

		free(inBuff);
	}

	printf("Starting relaxation for task %d.\n", taskID);
	// Do the work
	relaxation(currArray, nextArray, chunkSize, arraySize, precision);

	printf("Finished relaxation for task %d.\n", taskID);

	if (taskID != MASTER)
	{
		printf("Sending data to master (task %d).\n", taskID);
		double* outData = malloc(chunkSize*arraySize*sizeof(double));
		int j, currPos = 0;
		for( j = 2; j < chunkSize; j++) // skip first 2 columns
		{
			memcpy(&outData[currPos], nextArray[j], arraySize*sizeof(double));
			currPos += arraySize;
		}

		MPI_Send(outData, arraySize*chunkSize, MPI_DOUBLE, MASTER, TAG_COMPLETE_DATA, MPI_COMM_WORLD);
		free(outData);
	}
	else
	{
		printf("Receiving data.\n");
		int i, columnPos = chunkSize;
		double* inData = malloc(chunkSize*arraySize*sizeof(double));
		for( i = 1; i < numTasks; i++)
		{
			chunkSize = calculateChunkSize(arraySize, i);

			MPI_Status status;
			MPI_Recv(inData, arraySize*chunkSize, MPI_DOUBLE, i, TAG_COMPLETE_DATA, MPI_COMM_WORLD, &status);

			int currPos = 0, j;
			for( j = 0; j < chunkSize - 2; j++)
			{
				memcpy(nextArray[columnPos++], &inData[currPos], arraySize*sizeof(double));
				currPos += arraySize;
			}
		}
		free(inData);
	}


	if (__VERBOSE && (taskID == MASTER))
	{
		printf("Parallel result:\n");
		printSquareArray(nextArray, arraySize);
	}

	if (testRight && (taskID == MASTER))
	{
		double** startArray	= make2DDoubleArray(arraySize, arraySize);
		double** endArray	= make2DDoubleArray(arraySize, arraySize);
		initArrayPattern(startArray, arraySize);

		int fin = 0;

		while( fin == 0 )
		{
			fin = averageFour(startArray, endArray, arraySize, arraySize, precision);

			int i;
			for( i = 0; i < arraySize - 1; i++)
			{
				memcpy(startArray[i], endArray[i], arraySize * sizeof(double));
			}
		}

		int diff = checkDiff(endArray, nextArray, arraySize, arraySize, precision);

		printf("Parallel result %s the serial result\n", diff ? "matched" : "didn't match");
	}

	printf("Starting cleanup for task %d.\n", taskID);

	// Cleanup
	free2DDoubleArray(currArray, arraySize);
	free2DDoubleArray(nextArray, arraySize);

	double eTime = MPI_Wtime();

	printf("Relaxed %d square matrix in %f seconds (task: %d/%d)\n",
			arraySize, eTime, taskID, numTasks);


	MPI_Finalize();

	system("pause");

	return 0;
}
