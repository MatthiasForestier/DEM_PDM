#pragma once
#include <vector>

class Logger {
public:
    std::vector<float> objectiveHistory;
    std::vector<float> gradientNormHistory;
    //std::vector<float> conditionHistory; 

    void logStep(float objectiveValue, float gradientNorm) { // float condition
        objectiveHistory.push_back(objectiveValue);
        gradientNormHistory.push_back(gradientNorm);
        //conditionHistory.push_back(condition);
    }
    
    // Optionally, a clear function.
    void clear() {
        objectiveHistory.clear();
        gradientNormHistory.clear();
        //conditionHistory.clear();
    }
};

