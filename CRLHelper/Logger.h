#pragma once
#include <vector>

class Logger {
public:
    // Buffers to store the evolution of data.
    std::vector<float> objectiveHistory;
    std::vector<float> gradientNormHistory;
    
    // Add a new log entry.
    void logStep(float objectiveValue, float gradientNorm) {
        objectiveHistory.push_back(objectiveValue);
        gradientNormHistory.push_back(gradientNorm);
    }
    
    // Optionally, add a method to clear logs.
    void clear() {
        objectiveHistory.clear();
        gradientNormHistory.clear();
    }
};
