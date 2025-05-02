#include <memory>
#include "TemplateApp.h"        // TemplateApp will also include SplashScreenResult.h
#include "CRLHelper/VecMatDef.h"


// The existing function (if still needed)
SplashScreenResult showSplashScreen();

// New SplashScreen class that holds a TemplateApp instance.
class SplashScreen {
public:
    SplashScreen() {
        m_result = showSplashScreen();
        // Create the TemplateApp instance using the obtained result.
        m_templateApp = std::make_shared<TemplateApp>(m_result);
    }
    
    const SplashScreenResult& getResult() const { return m_result; }
    
    std::shared_ptr<TemplateApp> getTemplateApp() const { return m_templateApp; }
    
    // Forward TemplateApp calls.
    Optimization::OptimizationStatus energyMinimizationStepDyn() {
        return m_templateApp->energyMinimizationStepDyn();
    }
    
    Optimization::OptimizationStatus energyMinimizationStep() {
        return m_templateApp->energyMinimizationStep();
    }
    
private:
    SplashScreenResult m_result;
    std::shared_ptr<TemplateApp> m_templateApp;
};

