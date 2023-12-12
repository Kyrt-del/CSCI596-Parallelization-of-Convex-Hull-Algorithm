// This code generates n random points and writes them to input.txt
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>

using namespace std;

class Point
{
public:
    double x, y;
    Point(double x, double y)
    {
        this->x = x;
        this->y = y;
    }
};

int main(int argc, char *argv[])
{
    // generate n random points
    int n = atoi(argv[1]);
    vector<Point> points;
    srand(time(NULL));
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distribution(100000, 100000);
    for (int i = 0; i < n; i++)
    {
        double x = distribution(generator);
        double y = distribution(generator);
        points.push_back(Point(x, y));
    }

    // write the points to input.txt
    ofstream outputFile("input.txt");
    if (!outputFile.is_open())
    {
        std::cerr << "Error opening the file." << std::endl;
        return 1; // Return an error code
    }
    outputFile << n << endl;
    for (auto point : points)
    {
        outputFile << point.x << " " << point.y << endl;
    }
    outputFile.close();

    return 0;
}