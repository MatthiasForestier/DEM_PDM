#include "CRLHelper/CRLApp.h"
#include "Projects/Template/include/App/TemplateApp.h"

int main() {
    CRLApp app;
    app.subapps_selector.second.push_back(std::make_shared<TemplateApp>());
    app.launch();

    return 0;
}
