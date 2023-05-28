#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

const double START_GAS = 20.0;

struct Node {
    int number;
    double price;
    double x;
    double y;
};

std::vector<Node> graph = {
    {0, 0.0, 0.0, 0.0},  // Placeholder for start node
    {1, 0.0, 0.0, 0.0},  // Placeholder for finish node
};

double randomUniform(double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(min, max);
    return dis(gen);
}

void initializeGraph() {
    const int numNodes = 21;  // Number of nodes excluding start and finish

    for (int i = 0; i < numNodes; ++i) {
        double price = randomUniform(1.0, 3.0);
        double x = randomUniform(200.0, 1000.0);
        double y = randomUniform(200.0, 1000.0);
        graph.push_back({i, price, x, y});
    }

    // Placeholder for start and finish nodes with random coordinates
    double startX = randomUniform(200.0, 1000.0);
    double startY = randomUniform(200.0, 1000.0);
    double finishX = randomUniform(200.0, 1000.0);
    double finishY = randomUniform(200.0, 1000.0);
    graph[0] = {-1, 0.0, startX, startY};      // Start node
    graph[1] = {-2, NAN, finishX, finishY};    // Finish node
}

std::vector<std::vector<double>> calculateDistances() {
    const int numNodes = graph.size() - 2;  // Number of nodes excluding start and finish
    std::vector<std::vector<double>> distances(numNodes, std::vector<double>(numNodes, 0.0));

    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            double dx = graph[i + 2].x - graph[j + 2].x;
            double dy = graph[i + 2].y - graph[j + 2].y;
            distances[i][j] = std::sqrt(dx * dx + dy * dy);
        }
    }

    return distances;
}

std::vector<std::vector<double>> calculateRangeValues(const std::vector<std::vector<double>>& distances) {
    const double Kilo_per_gas = 10.0;
    const double START_GAS = 20.0;
    const double MAX_GAS = 40.0;
    const double START_RANGE = START_GAS * Kilo_per_gas;
    const double MAX_RANGE = MAX_GAS * Kilo_per_gas;
    const int numPoints = graph.size() - 2;  // Number of nodes excluding start and finish

    std::vector<std::vector<double>> rangeValues(numPoints);

    for (int u = 0; u < numPoints; ++u) {
        rangeValues[u].reserve(numPoints);
        double currentRange = MAX_RANGE;

        for (int v = 0; v < numPoints; ++v) {
            if (graph[v + 2].price < graph[u + 2].price && distances[u][v] <= MAX_RANGE) {
                double range = MAX_RANGE - distances[u][v];
                rangeValues[u].push_back(range);
            }
        }

        rangeValues[u].push_back(0.0);  // Add 0 as a possible range value
        std::sort(rangeValues[u].begin(), rangeValues[u].end());
    }

    return rangeValues;
}

std::vector<std::vector<std::pair<double, std::vector<int>>>> initializeCostDict(
    const std::vector<std::vector<double>>& distances,
    const std::vector<std::vector<double>>& rangeValues
) {
    const double Kilo_per_gas = 10.0;
    const double MAX_RANGE = 40.0 * Kilo_per_gas;
    const int numPoints = graph.size() - 2;  // Number of nodes excluding start and finish
    
    std::vector<std::vector<std::pair<double, std::vector<int>>>> cost(
        numPoints,
        std::vector<std::pair<double, std::vector<int>>>(rangeValues[0].size())
    );

    for (int u = 0; u < numPoints; ++u) {
        for (int rvIdx = 0; rvIdx < rangeValues[u].size(); ++rvIdx) {
            double rv = rangeValues[u][rvIdx];
            //std::cout<< rv<< std::endl;
            if (rv <= MAX_RANGE - distances[u][numPoints - 1]) {
                double c = (distances[u][numPoints - 1] - rv) * (graph[u + 2].price / Kilo_per_gas);
                cost[u][rvIdx] = {c, {numPoints}};
            } else {
                cost[u][rvIdx] = {INFINITY, {}};
            }
        }
    }

    return cost;
}

std::vector<std::vector<std::pair<double, std::vector<int>>>> calculateFinalCost(
    const std::vector<std::vector<double>>& distances,
    const std::vector<std::vector<double>>& rangeValues,
    const std::vector<std::vector<std::pair<double, std::vector<int>>>>& oldCost
) {
    const double Kilo_per_gas = 10.0;
    const double MAX_RANGE = 40.0 * Kilo_per_gas;
    const int numPoints = graph.size() - 2;  // Number of nodes excluding start and finish

    std::vector<std::vector<std::pair<double, std::vector<int>>>> cost(
        numPoints,
        std::vector<std::pair<double, std::vector<int>>>(rangeValues[0].size())
    );

    for (int u = 0; u < numPoints; ++u) {
        std::vector<std::pair<double, std::vector<int>>> indep1v;
        std::vector<std::pair<double, std::vector<int>>> indep2v;

        for (int v = 0; v < numPoints; ++v) {
            double d = distances[u][v];

            if (d > MAX_RANGE || u == v) {
                continue;
            }

            if (graph[v + 2].price <= graph[u + 2].price) {
                const auto& [c, prev] = oldCost[v][0];

                if (!prev.empty()) {
                    std::vector<int> newPrev = prev;
                    newPrev.push_back(v);
                    indep1v.push_back({c + d * (graph[u + 2].price / Kilo_per_gas), newPrev});
                }
            } else {
                const auto& [c, prev] = oldCost[v][rangeValues[v].size() - 1];

                if (!prev.empty()) {
                    std::vector<int> newPrev = prev;
                    newPrev.push_back(v);
                    indep2v.push_back({c + MAX_RANGE * (graph[u + 2].price / Kilo_per_gas), newPrev});
                }
            }
        }

        std::sort(indep1v.rbegin(), indep1v.rend());

        std::pair<double, std::vector<int>> min_indep2v;
        if (!indep2v.empty()) {
            min_indep2v = *std::min_element(indep2v.begin(), indep2v.end());
        } else {
            min_indep2v = {INFINITY, {}};
        }

        bool from_indep2v = true;
        if (!indep1v.empty()) {
            from_indep2v = min_indep2v <= indep1v.back();
        }

        for (int rvIdx = 0; rvIdx < rangeValues[u].size(); ++rvIdx) {
            double rv = rangeValues[u][rvIdx];

            if (from_indep2v) {
                cost[u][rvIdx] = {min_indep2v.first - rv * (graph[u + 2].price / Kilo_per_gas), min_indep2v.second};
            } else {
                while (!indep1v.empty() && distances[indep1v.back().second.back()][u] < rv) {
                    indep1v.pop_back();
                }

                if (!indep1v.empty() && min_indep2v > indep1v.back()) {
                    cost[u][rvIdx] = {indep1v.back().first - rv * (graph[u + 2].price / Kilo_per_gas), indep1v.back().second};
                } else {
                    from_indep2v = true;
                    cost[u][rvIdx] = {min_indep2v.first - rv * (graph[u + 2].price / Kilo_per_gas), min_indep2v.second};
                }
            }
        }
    }

    return cost;
}

std::pair<double, std::vector<int>> findStartFinishCost(
    const std::vector<std::vector<double>>& distances,
    const std::vector<std::vector<std::pair<double, std::vector<int>>>>& cost
) {
    const double Kilo_per_gas = 10.0;
    const double START_RANGE = 20.0 * Kilo_per_gas;
    const int numPoints = graph.size() - 2;  // Number of nodes excluding start and finish

    double minCost = INFINITY;
    std::vector<int> prev;
    for (int v = 0; v < numPoints; ++v) {
        double d = distances[v][numPoints];

        if (d > START_RANGE) {
            continue;
        }

        const auto& [cCost, prevV] = cost[v][cost[v].size() - 1 - (int)(START_RANGE - d)];
        if (cCost <= minCost) {
            minCost = cCost;
            prev = prevV;
            prev.push_back(v);
            prev.erase(prev.begin());  // Pop the finishing index
        }
    }

    return {minCost, prev};
}

void printPath( double minCost, const std::vector<int>& prev) {
    if (minCost != INFINITY) {
        std::cout << "Path found with cost: " << minCost << "!" << std::endl;
        std::cout << "start ==> ";
        for (int i = prev.size() - 1; i >= 0; --i) {
            std::cout << graph[prev[i] + 2].number << " ==> ";
        }
        std::cout << "finish" << std::endl;
    } else {
        std::cout << "No path found." << std::endl;
    }
}

int main() {
    initializeGraph();
    std::vector<std::vector<double>> distances = calculateDistances(); //working
    
    // for (auto i: distances)
    // {
    //     for (auto j: i)
    //     {
    //         std::cout<<j<< std::endl;
    //     }
    // }

    std::vector<std::vector<double>> rangeValues = calculateRangeValues(distances); //working
    // for (auto i: rangeValues)
    // {
    //     for (auto j: i)
    //     {
    //         std::cout<<j<< std::endl;
    //     }
    // }
    std::vector<std::vector<std::pair<double, std::vector<int>>>> cost = initializeCostDict(distances, rangeValues);
    
    for (auto i: cost)
    {
        for (auto j: i)
        {
            std::cout<<j.first<<" : ";
            for(auto k: j.second)
            {
                std::cout<<k<<std::endl;
            }

        }
    }
    // for (int i = 0; i < 20; ++i) {
    //     cost = calculateFinalCost(distances, rangeValues, cost);
    // }

    // auto [minCost, prev] = findStartFinishCost(distances, cost);
    // printPath(minCost, prev);

    return 0;
}
