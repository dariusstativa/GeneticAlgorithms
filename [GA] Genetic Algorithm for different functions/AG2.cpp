#include <iostream>
#include <vector>
#include <bitset>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <iomanip> 

using namespace std;

const int BITS_PER_DIM = 20;      
const int MAX_VALUE = (1 << BITS_PER_DIM) - 1;
double v[31],k=0;

const int POPULATION_SIZE = 100; 
const int GENERATIONS = 2000;     
const double MUTATION_RATE = 0.01;  
const double CROSSOVER_RATE = 0.5;  
const int ELITISM_COUNT = 5;       


double LOWER_BOUND, UPPER_BOUND;


double rastrigin(const vector<double>& x) {
    double result = 10.0 * x.size();
    for (double xi : x) {
        result += xi * xi - 10.0 * cos(2 * M_PI * xi);
    }
    return result;
}


double schwefel(const vector<double>& x) {
    double result = 0.0;
    for (double xi : x) {
        result += -xi * sin(sqrt(abs(xi)));
    }
    return result;
}


double michalewicz(const vector<double>& x, double m = 10.0) {
    double result = 0.0;
    for (int i = 0; i < x.size(); i++) {
        result -= sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / M_PI), 2 * m);
    }
    return result;
}


double deJong(const vector<double>& x) {
    double result = 0.0;
    for (double xi : x) {
        result += xi * xi;
    }
    return result;
}



double decode(const bitset<BITS_PER_DIM>& bits) {
    int integer_value = (int)bits.to_ulong();
    double ratio = (double)integer_value / MAX_VALUE;
    double decodedValue = LOWER_BOUND + ratio * (UPPER_BOUND - LOWER_BOUND);
    return decodedValue;
}


vector<bitset<BITS_PER_DIM>> generateRandomSolution(int dimensions) {
    vector<bitset<BITS_PER_DIM>> solution(dimensions);
    for (int i = 0; i < dimensions; i++) {
        solution[i] = bitset<BITS_PER_DIM>(rand() % (MAX_VALUE + 1));
    }
    return solution;
}


vector<double> decodeSolution(const vector<bitset<BITS_PER_DIM>>& solution) {
    vector<double> decoded(solution.size());
    for (size_t i = 0; i < solution.size(); i++) {
        decoded[i] = decode(solution[i]);
    }
    return decoded;
}

pair<vector<bitset<BITS_PER_DIM>>, vector<bitset<BITS_PER_DIM>>> uniformCrossover(
    const vector<bitset<BITS_PER_DIM>>& parent1,
    const vector<bitset<BITS_PER_DIM>>& parent2) {
    vector<bitset<BITS_PER_DIM>> offspring1 = parent1;
    vector<bitset<BITS_PER_DIM>> offspring2 = parent2;

    if (((double)rand() / RAND_MAX) > CROSSOVER_RATE) {
        for (size_t dim = 0; dim < parent1.size(); dim++) {
            for (int i = 0; i < BITS_PER_DIM; i++) {
                if (rand() % 2) {
                    offspring1[dim][i] = parent2[dim][i];
                    offspring2[dim][i] = parent1[dim][i];
                }
            }
        }
    }
    return {offspring1, offspring2};
}


void mutate(vector<bitset<BITS_PER_DIM>>& individual) {
    for (size_t dim = 0; dim < individual.size(); dim++) {
        for (int i = 0; i < BITS_PER_DIM; i++) {
            if (((double)rand() / RAND_MAX) < MUTATION_RATE) {
                individual[dim].flip(i);
            }
        }
    }
}


vector<vector<bitset<BITS_PER_DIM>>> tournamentSelection(const vector<vector<bitset<BITS_PER_DIM>>>& population,
                                                         const vector<double>& fitness) {
    vector<vector<bitset<BITS_PER_DIM>>> selected;
    while (selected.size() < POPULATION_SIZE) {
        vector<int> candidates(TOURNAMENT_SIZE);
        for (int i = 0; i < TOURNAMENT_SIZE; i++) {
            candidates[i] = rand() % POPULATION_SIZE;
        }

        int bestCandidate = candidates[0];
        for (int i = 1; i < TOURNAMENT_SIZE; i++) {
            if (fitness[candidates[i]] < fitness[bestCandidate]) {
                bestCandidate = candidates[i];
            }
        }
        selected.push_back(population[bestCandidate]);
    }
    return selected;
}


void geneticAlgorithm(int dimensions, int functionChoice) {
    srand(time(0));

   
    vector<vector<bitset<BITS_PER_DIM>>> population(POPULATION_SIZE);
    for (int i = 0; i < POPULATION_SIZE; i++) {
        population[i] = generateRandomSolution(dimensions);
    }

 
    vector<double> fitness(POPULATION_SIZE);
    for (int i = 0; i < POPULATION_SIZE; i++) {
        vector<double> decoded = decodeSolution(population[i]);
        switch (functionChoice) {
            case 1: fitness[i] = rastrigin(decoded); break;
            case 2: fitness[i] = schwefel(decoded); break;
            case 3: fitness[i] = michalewicz(decoded); break;
            case 4: fitness[i] = deJong(decoded); break;
        }
    }

    double bestFitness = numeric_limits<double>::max();
    vector<double> bestSolution;

    for (int gen = 0; gen < GENERATIONS; gen++) {
        vector<vector<bitset<BITS_PER_DIM>>> newPopulation;

      
        vector<int> eliteIndices(POPULATION_SIZE);
        iota(eliteIndices.begin(), eliteIndices.end(), 0);
        sort(eliteIndices.begin(), eliteIndices.end(), [&](int a, int b) {
            return fitness[a] < fitness[b];
        });

        for (int i = 0; i < ELITISM_COUNT; i++) {
            newPopulation.push_back(population[eliteIndices[i]]);
        }

       
        while (newPopulation.size() < POPULATION_SIZE) {
            vector<bitset<BITS_PER_DIM>> parent1 = tournamentSelection(population, fitness)[0];
            vector<bitset<BITS_PER_DIM>> parent2 = tournamentSelection(population, fitness)[0];

            auto [offspring1, offspring2] = uniformCrossover(parent1, parent2);

            mutate(offspring1);
            mutate(offspring2);

            newPopulation.push_back(offspring1);
            if (newPopulation.size() < POPULATION_SIZE) {
                newPopulation.push_back(offspring2);
            }
        }

        population = newPopulation;

      
        for (int i = 0; i < POPULATION_SIZE; i++) {
            vector<double> decoded = decodeSolution(population[i]);
            switch (functionChoice) {
                case 1: fitness[i] = rastrigin(decoded); break;
                case 2: fitness[i] = schwefel(decoded); break;
                case 3: fitness[i] = michalewicz(decoded); break;
                case 4: fitness[i] = deJong(decoded); break;
            }
        }

      
        auto currentBestIt = min_element(fitness.begin(), fitness.end());
        double currentBestFitness = *currentBestIt;
        if (currentBestFitness < bestFitness) {
            bestFitness = currentBestFitness;
            bestSolution = decodeSolution(population[currentBestIt - fitness.begin()]);
        }

        
    }

    cout << "\nBest fitness found: " << fixed << setprecision(5) << bestFitness << endl;
    
    cout << "Solution: ";
    for (double val : bestSolution) {
        cout << fixed << setprecision(5) << val << " ";
    }
    cout << endl;
}

int main() {
    int dimensions, functionChoice;

    cout << "Enter the problem dimensions (5, 10, 30): ";
    cin >> dimensions;

    cout << "Choose the function to optimize:\n";
    cout << "1. Rastrigin\n";
    cout << "2. Schwefel\n";
    cout << "3. Michalewicz\n";
    cout << "4. De Jong (Sum of Squares)\n";
    cin >> functionChoice;

    switch (functionChoice) {
        case 1: LOWER_BOUND = -5.12; UPPER_BOUND = 5.12; break;
        case 2: LOWER_BOUND = -500; UPPER_BOUND = 500; break;
        case 3: LOWER_BOUND = 0; UPPER_BOUND = M_PI; break;
        case 4: LOWER_BOUND = -5.12; UPPER_BOUND = 5.12; break;
        default: cout << "Invalid choice! Exiting.\n"; return 0;
    }

    
for(int i=1;i<=30;i++)
geneticAlgorithm(dimensions, functionChoice);
    return 0;
}
