// SplashScreen.cpp
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "Projects/Template/include/App/SplashScreen.h"

SplashScreenResult showSplashScreen() {
    SplashScreenResult result;
    // Default values for interactive mode:
    result.use3D = false;
    result.exportMode = false;
    // Default values for export mode parameters:
    result.timeStep = 0.001;
    result.dynamic = true;
    result.viscosity = false;
    result.numSimulations = 5;
    result.simulationTime = 3.0;
    result.animationStartTime = 6.0;
    result.exponent_convergence_threshold = -3;
    result.maxIter = 3000;
    result.numParticles = 40;
    result.V0 = 0.5;
    result.L = 0.5;
    result.densify = true;
    bool done = false;

    if (!glfwInit()) {
        return result;
    }
    
    GLFWwindow* window = glfwCreateWindow(1400, 600, "Select Mode", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return result;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;
    io.FontGlobalScale = 3;

    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    while (!glfwWindowShouldClose(window) && !done) {
        glfwPollEvents();
    
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
    
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
    
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2((float)display_w, (float)display_h));
    
        ImGui::Begin("Select Mode", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoTitleBar);
        ImGui::Text("Choose application mode:");

        // Mode selection: interactive or export
        static int mode = 0; // 0 = Interactive, 1 = Export Simulation
        ImGui::RadioButton("Interactive Mode", &mode, 0);
        ImGui::RadioButton("Export Simulation", &mode, 1);
        
        if (mode == 0) {
            result.exportMode = false;
            // Show the two classic buttons for 2D and 3D:
            float availWidth = ImGui::GetContentRegionAvail().x;
            float buttonWidth = (availWidth - ImGui::GetStyle().ItemSpacing.x) * 0.5f;
            if (ImGui::Button("2D Mode", ImVec2(buttonWidth, 0))) {
                result.use3D = false;
                done = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("3D Mode", ImVec2(buttonWidth, 0))) {
                result.use3D = true;
                done = true;
            }
        } else {
            result.exportMode = true;
            // Inside export mode branch:
            ImGui::InputDouble("Time Step", &result.timeStep, 0.001, 0.01, "%.3f");
            ImGui::Checkbox("Dynamic", &result.dynamic);
            ImGui::InputInt("Epsilon Dynamic", &result.exponent_convergence_threshold);
            ImGui::InputDouble("Tolerable iterations", &result.maxIter);
            ImGui::Checkbox("Viscosity", &result.viscosity);
            ImGui::InputInt("Number of Simulations", &result.numSimulations);
            ImGui::InputInt("Number of Particles", &result.numParticles);
            ImGui::Checkbox("Increase density", &result.densify);
            ImGui::InputDouble("Simulation Duration (s)", &result.simulationTime, 1.0, 5.0, "%.1f");
            ImGui::InputDouble("Animation Start (s)", &result.animationStartTime, 0.1, 1.0, "%.2f");

            // Add scenario shape selection. Circle is default (index 0).
            static int shapeIndex = 0;
            const char* scenarioShapes[] = {"Circle","Square","Tunnel"};
            ImGui::Combo("Scenario Object",&shapeIndex,scenarioShapes,IM_ARRAYSIZE(scenarioShapes));
            result.scenarioShapeIndex = shapeIndex;
        
            // only for tunnel do we need shear / periodic parameters
            if (shapeIndex == 2) {
                ImGui::Checkbox   ("Periodic X",       &result.periodicX);
                ImGui::InputDouble("V0 (max speed)",    &result.V0,           0.1, 1.0, "%.2f");
                ImGui::InputDouble("L (half-height)",   &result.L,            0.1, 1.0, "%.2f");
            }
            result.scenarioShapeIndex = shapeIndex;

            if (ImGui::Button("Run Export Simulation")) {
                done = true;
            }
        }
        ImGui::End();
    
        ImGui::Render();
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }
    
    // Clean up splash screen resources.
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    
    return result;
}