// #include <H5Cpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>


namespace stdfs = std::filesystem;

#include "Projects/Template/include/App/SplashScreen.h"
#include "Projects/Template/include/Model/MassSpring.h"
#include "Projects/Template/include/App/TemplateApp.h"


// Runs the simulation in export mode and writes behavior to a hard-coded text file.
void exportSimulationForML(const SplashScreenResult &params);

//void exportSimulationForMLHD5(const SplashScreenResult &params);
