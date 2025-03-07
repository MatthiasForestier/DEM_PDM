#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Projects/Template/include/App/SplashScreen.h"

// Function that shows a splash screen to let the user choose between 2D and 3D.
// Returns true for 3D mode, false for 2D mode.
SplashScreenResult showSplashScreen() {
    SplashScreenResult result;
    result.use3D = true;        // default mode: 3D
    // result.numParticles = 100;  // default number of particles

    bool done = false;

    // Initialize GLFW
    if (!glfwInit()) {
        return result;
    }

    // Create a small window for the splash screen
    GLFWwindow* window = glfwCreateWindow(800, 400, "Select Mode", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return result;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);  // Enable vsync

    // Setup ImGui context and initialize for GLFW and OpenGL3
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;
    io.FontGlobalScale = 3;  // Increase font size

    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    // Splash screen loop: run until the window is closed or the user makes a selection.
    while (!glfwWindowShouldClose(window) && !done) {
        glfwPollEvents();
    
        // Start new ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
    
        // Get the current framebuffer size from GLFW
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
    
        // Set the next window's position and size to match the entire GLFW window
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2((float)display_w, (float)display_h));
    
        // Create an ImGui window without a title bar and without resize
        ImGui::Begin("Select Mode", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoTitleBar);
        ImGui::Text("Choose application mode:");
    
        // Let the user input the number of particles
        // ImGui::InputInt("Number of Particles", &result.numParticles);
        // // Ensure a valid (positive) number is selected
        // if (result.numParticles < 1) {
        //     result.numParticles = 1;
        // }
    
        // Compute available width for buttons and create two buttons that fill the width equally.
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
        ImGui::End();
    
        // Render the splash screen
        ImGui::Render();
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }
    
    // Clean up resources for the splash screen
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    return result;
}
