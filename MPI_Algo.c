#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// Function to determine the orientation of three points (p, q, r)
int orientation(const double *p, const double *q, const double *r)
{
    double val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);
    if (val == 0.0)
        return 0;             // Collinear
    return (val > 0) ? 1 : 2; // Clockwise or Counterclockwise
}

// Function comparePoints is used to sort the polar angles array
int comparePoints(const void *a, const void *b)
{
    const double *p1 = *(const double **)a;
    const double *p2 = *(const double **)b;
    return (p1[0] < p2[0]) || (p1[0] == p2[0] && p1[1] < p2[1]) ? -1 : 1;
}

// Function to compute the Convex Hull using Graham's Scan algorithm
void computeConvexHull(const double *points, int size, double **convexHull, int *convexHullSize)
{
    // Find the leftmost point in points
    int leftmost = 0;
    int leftmost_count = 1;
    for (int i = 1; i < size; i++)
    {
        if (points[i * 2] < points[leftmost * 2])
        {
            leftmost = i;
            leftmost_count = 1;
        }
        else if (points[i * 2] == points[leftmost * 2] && points[i * 2 + 1] < points[leftmost * 2 + 1])
        {
            leftmost_count += 1;
        }
    }

    // make an array of array containing polar angles of all points from leftmost point
    int polarAnglesSize = size - leftmost_count;
    double **polarAngles = malloc(polarAnglesSize * sizeof(double *));
    int rem = 0;
    for (int i = 0; i < size; i++)
    {
        if (points[i * 2] == points[leftmost * 2] && points[i * 2 + 1] == points[leftmost * 2 + 1])
        {
            rem++;
            continue;
        }
        polarAngles[i - rem] = malloc(2 * sizeof(double));
        polarAngles[i - rem][0] = atan((points[i * 2 + 1] - points[leftmost * 2 + 1]) / (points[i * 2] - points[leftmost * 2]));
        // polarAngles[i - rem][0] = atan2(points[i * 2 + 1] - points[leftmost * 2 + 1], points[i * 2] - points[leftmost * 2]);
        polarAngles[i - rem][1] = i;
    }

    // sort the polar angles array based on polar angles
    qsort(polarAngles, polarAnglesSize, sizeof(double *), comparePoints);

    // make a stack and push the leftmost point and the point with the smallest polar angle
    int *stack = malloc(size * sizeof(int));
    int top = -1;
    stack[++top] = leftmost;
    stack[++top] = (int)polarAngles[0][1];
    // compute the Convex Hull
    for (int i = 1; i < polarAnglesSize; i++)
    {
        int p1 = stack[top - 1];
        int p2 = stack[top];
        int p3 = polarAngles[i][1];
        double crossProduct = (points[p2 * 2] - points[p1 * 2]) * (points[p3 * 2 + 1] - points[p1 * 2 + 1]) - (points[p2 * 2 + 1] - points[p1 * 2 + 1]) * (points[p3 * 2] - points[p1 * 2]);
        while (crossProduct <= 0)
        {
            top--;
            p1 = stack[top - 1];
            p2 = stack[top];
            crossProduct = (points[p2 * 2] - points[p1 * 2]) * (points[p3 * 2 + 1] - points[p1 * 2 + 1]) - (points[p2 * 2 + 1] - points[p1 * 2 + 1]) * (points[p3 * 2] - points[p1 * 2]);
        }
        stack[++top] = p3;
    }

    // Copy the Convex Hull to the output array
    *convexHullSize = top + 1;
    *convexHull = malloc((*convexHullSize) * 2 * sizeof(double));
    for (int i = 0; i <= top; i++)
    {
        (*convexHull)[i * 2] = points[(int)stack[i] * 2];
        (*convexHull)[i * 2 + 1] = points[(int)stack[i] * 2 + 1];
    }

    free(stack);
    free(polarAngles);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int totalPoints = 0;
    double *allPoints = NULL;

    if (rank == 0)
    {
        //  Master process reads points from file
        FILE *file = fopen("input.txt", "r");
        if (!file)
        {
            fprintf(stderr, "Error opening the file.\n");
            // MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Count the total number of points
        fscanf(file, "%d", &totalPoints);

        // Allocate memory for all points
        allPoints = malloc(totalPoints * 2 * sizeof(double));

        // Read points from file
        for (int i = 0; i < totalPoints; i++)
        {
            fscanf(file, "%lf %lf", &allPoints[i * 2], &allPoints[i * 2 + 1]);
        }

        fclose(file);
    }

    // Broadcast the total number of points to all processes
    MPI_Bcast(&totalPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Scatter the array of points to all MPI processes
    int pointsPerProcess = totalPoints / size;
    double *localPoints = malloc(pointsPerProcess * 2 * sizeof(double));
    MPI_Scatter(allPoints, pointsPerProcess * 2, MPI_DOUBLE, localPoints, pointsPerProcess * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Each MPI process computes the Convex Hull for its subset of points
    double *localConvexHull;
    int localConvexHullSize;
    computeConvexHull(localPoints, pointsPerProcess, &localConvexHull, &localConvexHullSize);

    int *localConvexHullSizes = NULL;
    if (rank == 0)
    {
        localConvexHullSizes = malloc(size * sizeof(int));
    }

    // Gather the size of local convex hulls to the root process
    MPI_Gather(&localConvexHullSize, 1, MPI_INT,
               localConvexHullSizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Gather the contents of local convex hulls to the root process
    double *gatheredConvexHulls = NULL;
    int *displacements = NULL;
    int totalGatherSize;
    if (rank == 0)
    {
        // Calculate displacements for MPI_Gatherv
        displacements = malloc(size * sizeof(int));
        displacements[0] = 0;
        for (int i = 1; i < size; i++)
        {
            displacements[i] = displacements[i - 1] + localConvexHullSizes[i - 1] * 2;
        }

        // Calculate the total size of the gathered convex hulls
        totalGatherSize = displacements[size - 1] + localConvexHullSizes[size - 1] * 2;
        gatheredConvexHulls = malloc(totalGatherSize * sizeof(double));
        // multiply all local convex hull sizes by 2 in the localConvexHullSizes array
        for (int i = 0; i < size; i++)
        {
            localConvexHullSizes[i] *= 2;
        }
    }
    
    MPI_Gatherv(localConvexHull, localConvexHullSize * 2, MPI_DOUBLE,
                gatheredConvexHulls, localConvexHullSizes, displacements, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // Merge the convex hulls on the root process
    double *mergedConvexHull;
    int mergedConvexHullSize;
    if (rank == 0)
    {
        // printf("Merging convex hulls\n");
        // printf("Total points: %d\n", totalPoints);
        computeConvexHull(gatheredConvexHulls, totalGatherSize / 2, &mergedConvexHull, &mergedConvexHullSize);

        // Write the merged convex hull to the output file
        FILE *outputFile = fopen("output_merged.txt", "w");
        for (int i = 0; i < mergedConvexHullSize; i++)
        {
            fprintf(outputFile, "%lf %lf\n", mergedConvexHull[i * 2], mergedConvexHull[i * 2 + 1]);
        }
        fclose(outputFile);

        free(gatheredConvexHulls);
        free(localConvexHullSizes);
        free(displacements);
    }

    // Free allocated memory
    free(allPoints);
    free(localPoints);
    free(localConvexHull);

    MPI_Finalize();
    return 0;
}
