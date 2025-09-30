#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

const int A = 10;
const double lowerBoundRastrigin = -5.12;
const double upperBoundRastrigin = 5.12;
const double lowerBoundMich = 0;
const double upperBoundMich = M_PI;
const double lowerBoundSphere = -5.12;
const double upperBoundSphere = 5.12;
const double lowerBoundSumSquare = -10;
const double upperBoundSumSquare = 10;
const double stepSize = 0.1;
int maxIterations = 100000;

double rastrigin(const std::vector<double>& x, const int n) {
    double sum = A * n;
    for (int i = 0; i < n; i++) {
        sum += x[i] * x[i] - A * cos(2 * M_PI * x[i]);
    }
    return sum;
}

double michalewicz(const std::vector<double>& x, const int n) {
    double sum = 0;
    const int m = 10;
    for (int i = 0; i < n; i++) {
        sum += sin(x[i]) * pow(sin(((i + 1) * x[i] * x[i]) / M_PI), 2 * m);
    }
    return -sum;
}

double sphere(const std::vector<double>& x, const int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}

double sumSquare(const std::vector<double>& x, const int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += (i + 1) * x[i] * x[i];
    }
    return sum;
}

std::vector<double> generateNeighbor(const std::vector<double>& current, const int n, const double lowerBound, const double upperBound) {
    std::vector<double> neighbor = current;
    int changedIndex = rand() % n;
    double step = ((double)rand() / RAND_MAX) * 2 * stepSize - stepSize;
    neighbor[changedIndex] += step;
    if (neighbor[changedIndex] < lowerBound) neighbor[changedIndex] = lowerBound;
    if (neighbor[changedIndex] > upperBound) neighbor[changedIndex] = upperBound;
    return neighbor;
}

std::vector<double> hillClimbingBestImprovement(const std::vector<double>& initial, const int n, double (*func)(const std::vector<double>&, const int), const double lowerBound, const double upperBound) {
    std::vector<double> current = initial;
    double currentFitness = func(initial, n);

    for (int i = 0; i < maxIterations; i++) {
        std::vector<double> neighbor = generateNeighbor(current, n, lowerBound, upperBound);
        double neighborFitness = func(neighbor, n);
        if (neighborFitness < currentFitness) {
            current = neighbor;
            currentFitness = neighborFitness;
        }
    }
    return current;
}

std::vector<double> hillClimbingFirstImprovement(const std::vector<double>& initial, const int n, double (*func)(const std::vector<double>&, const int), const double lowerBound, const double upperBound) {
    std::vector<double> current = initial;
    double currentFitness = func(initial, n);

    for (int i = 0; i < maxIterations; i++) {
        std::vector<double> neighbor = generateNeighbor(current, n, lowerBound, upperBound);
        double neighborFitness = func(neighbor, n);
        if (neighborFitness < currentFitness) {
            current = neighbor;
            currentFitness = neighborFitness;
        }
    }
    return current;
}

std::vector<double> hillClimbingWorstImprovement(const std::vector<double>& initial, const int n, double (*func)(const std::vector<double>&, const int), const double lowerBound, const double upperBound) {
    std::vector<double> current = initial;
    double currentFitness = func(initial, n);

    for (int i = 0; i < maxIterations; i++) {
        std::vector<double> worstNeighbor = current;
        double worstNeighborFitness = currentFitness;

        for (int j = 0; j < n; j++) {
            std::vector<double> neighbor = generateNeighbor(current, n, lowerBound, upperBound);
            double neighborFitness = func(neighbor, n);
            if (neighborFitness > worstNeighborFitness && neighborFitness < currentFitness) {
                worstNeighbor = neighbor;
                worstNeighborFitness = neighborFitness;
            }
        }

        if (worstNeighborFitness < currentFitness) {
            current = worstNeighbor;
            currentFitness = worstNeighborFitness;
        } else {
            break;
        }
    }
    return current;
}

std::vector<double> hillClimbingSimulatedAnnealing(const std::vector<double>& initial, const int n, double (*func)(const std::vector<double>&, const int), const double lowerBound, const double upperBound) {
    std::vector<double> current = initial;
    double currentFitness = func(initial, n);
    double temperature = 1000.0;
    double coolingRate = 0.003;

    for (int i = 0; i < maxIterations; i++) {
        std::vector<double> neighbor = generateNeighbor(current, n, lowerBound, upperBound);
        double neighborFitness = func(neighbor, n);

        if (neighborFitness < currentFitness) {
            current = neighbor;
            currentFitness = neighborFitness;
        } else {
            double acceptanceProbability = exp((currentFitness - neighborFitness) / temperature);
            if (((double)rand() / RAND_MAX) < acceptanceProbability) {
                current = neighbor;
                currentFitness = neighborFitness;
            }
        }

        temperature *= 1 - coolingRate;
        if (temperature < 1) break;
    }
    return current;
}

void testFunction(double (*func)(const std::vector<double>&, const int), double lowerBound, double upperBound, int n, const std::string& method) {
    double finalMinimum;
    std::vector<double> finalMinimumValues;
    std::vector<double> allResults;

    for (int j = 0; j < 100; j++) {
        std::vector<double> start(n);
        for (double& x : start) {
            x = lowerBound + ((double)rand() / RAND_MAX) * (upperBound - lowerBound);
        }

        std::vector<double> result;
        if (method == "Best Improvement") {
            result = hillClimbingBestImprovement(start, n, func, lowerBound, upperBound);
        } else if (method == "First Improvement") {
            result = hillClimbingFirstImprovement(start, n, func, lowerBound, upperBound);
        } else if (method == "Worst Improvement") {
            result = hillClimbingWorstImprovement(start, n, func, lowerBound, upperBound);
        } else if (method == "Simulated Annealing Hybrid") {
            result = hillClimbingSimulatedAnnealing(start, n, func, lowerBound, upperBound);
        }

        double resultValue = func(result, n);
        allResults.push_back(resultValue);

        if (j == 0 || resultValue < finalMinimum) {
            finalMinimum = resultValue;
            finalMinimumValues = result;
        }

        cout << "Best solution for current iteration (" << method << "): ";
        for (double x : result) {
            cout << x << " ";
        }
        cout << " -> Function value: " << resultValue << endl;
    }

    double mean = accumulate(allResults.begin(), allResults.end(), 0.0) / allResults.size();
    double stddev = sqrt(accumulate(allResults.begin(), allResults.end(), 0.0,
                     [mean](double acc, double val) { return acc + (val - mean) * (val - mean); }) / allResults.size());
    double minResult = *min_element(allResults.begin(), allResults.end());
    double maxResult = *max_element(allResults.begin(), allResults.end());

    cout << "\nFinal Best solution for " << method << ": ";
    for (double x : finalMinimumValues) {
        cout << x << " ";
    }
    cout << "\nFinal minimum function value: " << finalMinimum << endl;
    cout << "Mean of values: " << mean << endl;
    cout << "Standard Deviation of values: " << stddev << endl;
    cout << "Minimum value found: " << minResult << endl;
    cout << "Maximum value found: " << maxResult << "\n" << endl;
}

int main() {
    srand(time(0));
    int n = 5;
    int choice;
    string method;

    cout << "Choose the function to optimize:\n";
    cout << "1. Rastrigin\n";
    cout << "2. Michalewicz\n";
    cout << "3. Sphere\n";
    cout << "4. Sum Square\n";
    cin >> choice;

    cout << "Choose the method:\n";
    cout << "1. Best Improvement\n";
    cout << "2. First Improvement\n";
    cout << "3. Worst Improvement\n";
    cout << "4. Simulated Annealing Hybrid\n";
    int methodChoice;
    cin >> methodChoice;

    if (methodChoice == 1) {
        method = "Best Improvement";
    } else if (methodChoice == 2) {
        method = "First Improvement";
    } else if (methodChoice == 3) {
        method = "Worst Improvement";
    } else if (methodChoice == 4) {
        method = "Simulated Annealing Hybrid";
    } else {
        cout << "Invalid method choice.\n";
        return 1;
    }

    switch (choice) {
        case 1:
            testFunction(rastrigin, lowerBoundRastrigin, upperBoundRastrigin, n, method);
            break;
        case 2:
            testFunction(michalewicz, lowerBoundMich, upperBoundMich, n, method);
            break;
        case 3:
            testFunction(sphere, lowerBoundSphere, upperBoundSphere, n, method);
            break;
        case 4:
            testFunction(sumSquare, lowerBoundSumSquare, upperBoundSumSquare, n, method);
            break;
        default:
            cout << "Invalid function choice.\n";
            return 1;
    }

    return 0;
}

