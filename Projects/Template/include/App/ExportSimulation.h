// ExportSimulation.h
#ifndef EXPORT_SIMULATION_H
#define EXPORT_SIMULATION_H

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
#include "Projects/Template/include/App/ExportSimulation.h"

namespace stdfs = std::filesystem;

#include "SplashScreen.h"
#include "Projects/Template/include/Model/MassSpring.h"
#include "Projects/Template/include/App/ExportSimulation.h"
#include "Projects/Template/include/App/TemplateApp.h"


// Runs the simulation in export mode and writes behavior to a hard-coded text file.
void exportSimulationForML(const SplashScreenResult &params);

//void exportSimulationForMLHD5(const SplashScreenResult &params);

#endif