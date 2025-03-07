#include "CRLHelper/CRLApp.h"
#include "Projects/Template/include/App/TemplateApp.h"
#include "Projects/Template/include/App/SplashScreen.h"

int main() {
    SplashScreenResult InitParam = showSplashScreen();
    CRLApp app;
    auto myTemplateApp = std::make_shared<TemplateApp>(InitParam);
    app.subapps_selector.second.push_back(myTemplateApp);
    app.launch();

    return 0;
}
