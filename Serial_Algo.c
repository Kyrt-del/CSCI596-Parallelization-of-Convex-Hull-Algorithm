#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to determine the orientation of three points (p, q, r)
int orientation(const double *p, const double *q, const double *r)
{
    double val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);
    if (val == 0.0)
        return 0;             // Collinear
    return (val > 0) ? 1 : 2; // Clockwise or Counterclockwise
}

// Function to compare points for sorting
int comparePoints(const void *a, const void *b)
{
    const double *p1 = *(const double **)a;
    const double *p2 = *(const double **)b;
    return (p1[0] < p2[0]) || (p1[0] == p2[0] && p1[1] < p2[1]) ? -1 : 1;
}

// Function comparePoints2 is used to sort the polar angles array
int comparePoints2(const void *a, const void *b)
{
    const double *p1 = *(const double **)a;
    const double *p2 = *(const double **)b;
    return (p1[0] < p2[0]) || (p1[0] == p2[0] && p1[1] < p2[1]) ? -1 : 1;
}

// Function to compute the Convex Hull using Graham's Scan algorithm
void computeConvexHull2(const double *points, int size, double **convexHull, int *convexHullSize)
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
    qsort(polarAngles, polarAnglesSize, sizeof(double *), comparePoints2);

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

// Function to compute the Convex Hull using Graham's Scan algorithm
void computeConvexHull(const double *points, int size, double **convexHull, int *convexHullSize)
{

    if (size < 3)
    {
        // Convex Hull is not possible
        *convexHullSize = size;
        *convexHull = malloc(size * 2 * sizeof(double));
        for (int i = 0; i < size; i++)
        {
            (*convexHull)[i * 2] = points[i * 2];
            (*convexHull)[i * 2 + 1] = points[i * 2 + 1];
        }
        return;
    }
    // Sort points based on x-coordinates
    double **sortedPoints = malloc(size * sizeof(double *));
    for (int i = 0; i < size; i++)
    {
        sortedPoints[i] = malloc(2 * sizeof(double));
        sortedPoints[i][0] = points[i * 2];
        sortedPoints[i][1] = points[i * 2 + 1];
    }
    qsort(sortedPoints, size, sizeof(double *), comparePoints);
    // Initialize the Convex Hull
    double **stack = malloc(size * sizeof(double *));
    int top = -1;
    stack[++top] = sortedPoints[0];
    stack[++top] = sortedPoints[1];
    // Compute the Convex Hull
    for (int i = 2; i < size; i++)
    {
        while (top > 0 && orientation(stack[top - 1], stack[top], sortedPoints[i]) != 2)
        {
            top--;
        }
        stack[++top] = sortedPoints[i];
    }
    // Copy the Convex Hull to the output array
    *convexHullSize = top + 1;
    *convexHull = malloc((*convexHullSize) * 2 * sizeof(double));
    for (int i = 0; i <= top; i++)
    {
        (*convexHull)[i * 2] = stack[i][0];
        (*convexHull)[i * 2 + 1] = stack[i][1];
        free(stack[i]);
    }
    free(sortedPoints);

    free(stack);
}

int main(int argc, char **argv)
{
    int totalPoints = 0;
    double *allPoints = NULL;

    FILE *file = fopen("input.txt", "r");
    if (!file)
    {
        fprintf(stderr, "Error opening the file.\n");
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

    double *localConvexHull;
    int localConvexHullSize;

    computeConvexHull2(allPoints, totalPoints, &localConvexHull, &localConvexHullSize);

    // Write local convex hull points to output file
    char filename[20];
    // sprintf(filename, "output_%d.txt", rank);
    sprintf(filename, "output_serial.txt");
    FILE *outputFile = fopen(filename, "w");
    for (int i = 0; i < localConvexHullSize; i++)
    {
        fprintf(outputFile, "%lf %lf\n", localConvexHull[i * 2], localConvexHull[i * 2 + 1]);
    }
    fclose(outputFile);

    return 0;
}
