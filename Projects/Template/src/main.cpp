// main.cpp
#include "CRLHelper/CRLApp.h"
#include "./Projects/Template/include/App/SplashScreen.h"
#include "./Projects/Template/include/App/ExportSimulation.h"

int main() {
    SplashScreen splashScreen;
    SplashScreenResult initParam = splashScreen.getResult();

    if (initParam.exportMode) {
        // Run export simulation and then exit.
        exportSimulationForML(initParam);
        return 0;
    }

    // Otherwise, launch the interactive application.
    CRLApp app;
    // Retrieve the TemplateApp from the splash screen.
    app.subapps_selector.second.push_back(splashScreen.getTemplateApp());
    app.launch();

    return 0;
}