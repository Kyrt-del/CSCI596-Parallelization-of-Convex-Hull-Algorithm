# CSCI596: Parallelization of Convex Hull Algorithm

Convex hull finding algorithms involves finding a set of points from a given pool of points, such that
the set of points form the convex shape which holds the pool of points. The points found are also part of
the given pool of points. Convex hull finding algorithm can be used in a multitude of fields like computer
graphics, image processing, robotics, etc.

In this repository we have proposed a parallel algorithm using MPI which can be used to find convex hull
in a 2-D space. The corresponding serial algorithm is also given. The procedure used to prepare the parallel
algorithm has been discussed below, followed by the results. For checking whether the algorithm is working correctly,
we have prepared a visualizer which makes a scatter plot of all the points in the data and makes a convex hull to
check whether the prepared hull encloses the given set of points.

## Files

- `Serial_Algo.c` : Serial algorithm for finding Convex Hull
- `MPI_Algo.c`: Parallelized algorith for finding Convex Hull (Parallelized using MPI)
- `generate_random_points.cpp`: Code for generating random set of points which are used as dataset
  - The first line is the number of points in the file
  - The name of the file containing points is "input.txt"
  - The number of points for the dataset can be passed as an argument while running the code. For example, if we want
  the file to contain 1000 points then we can run the executable in the following way (assuming that the executable is named
  `generate_random_points`)

  ```sh
  ./generate_random_points 1000
  ```

- `visualizer.py`: Creates a scatter plot of points in the dataset and makes a convex hull using output file to check
if the convex hull found is correct.

## Implementation Details

### Serial Algorithm

We have used the Graham Scan algorithm for finding the convex hull in the serial algorithm. The steps in algorithm are discussed briefly
in the steps below

1. Find the leftmost point in terms of x-coordinate in the set of points.
2. Find the polar angles of all points with respect to leftmost point and sort the set of points (excluding the leftmost point) based on the polar angles found
3. Prepare a array to represent a stack and store the leftmost point and 1st point in the sorted array of points in the stack.
4. Start scan from the 2nd point in the sorted array.
    - Set p1 and p2 as the top points in the stack, and p3 as the current scanning point (p3 comes from the sorted array).
    - Assuming a straight line from p1 to p2, find if p3 lies of the left side or the right side by using the cross product.
    - If p3 is on the left then add it in the stack
    - If p1, p2 and p3 are collinear then remove p2 from the stack and start the scan again.
    - If p3 is on the right, then remove p2 from the stack and start the scan again.
5. Once the scan is finished, the points remaining in the stack are the points for convex hull.

### Parallel Algorithm

The natural observation while running serial algorithm for finding convex hulls is that the number of points which are used to define the
convex hull are far less than the total number of points in the dataset. And those points can be considered as the representation of the
whole dataset. We used this concept while formulating our methodology for parallelizing the algorithm. The steps for parallelization are:

1. Distribute the points in the dataset among all the threads equally using MPI_Scatter.
2. Each thread computes a convex hulls for the points sent to them.
3. All threads send the size of the convex hull or the number of points in the convex hulls to the master thread using MPI_Gather.
4. Once the sizes of local convex hulls are known to the master thread the displacement array is created in the master thread.
5. Gather all the local convex hulls in one array in master thread (called allConvexHulls) using the displacement array with the help of MPI_Gatherv call.
6. Find the convex hull using the allConvexHulls array in the master thread.

The process has been explained with the help of a diagram below for the case of 4 threads. The master thread's processes are colored green. All processes done by a thread are marked in the same color.

![mpi-algo-image]

## Results

The time vs thread curve is shown below.

![time-vs-thread-image]

The speedup vs thread curve is shown below.

![speedup-vs-thread-image]

The efficiency vs thread curve is shown below.

![efficiency-vs-thread-image]

`visualizer.py` can be used to visualize the convex hull prepared along with the input points. Below is an image of convex hull computed from
the dataset consisting of 10 million points using the proposed parallel algorithm. We can see that the convex hull is computed correctly.

![convex-hull-plot]

[mpi-algo-image]: images/Convex_Hull_Parallel.png
[time-vs-thread-image]: images/time_vs_threads.png
[speedup-vs-thread-image]: images/speedup_vs_threads.png
[efficiency-vs-thread-image]: images/efficiency_vs_threads.png
[convex-hull-plot]: images/Convex_hull_plot.png
