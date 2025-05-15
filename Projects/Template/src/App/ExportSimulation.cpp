#include "Projects/Template/include/App/ExportSimulation.h"

namespace stdfs = std::filesystem;

// Helper function to get the current date/time as a string in the format YYYY-MM-DD_HH-MM-SS.
std::string getCurrentDateTimeString() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    
    std::tm tm;
#ifdef _WIN32
    localtime_s(&tm, &now_time);
#else
    localtime_r(&now_time, &tm);
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    return oss.str();
}

void exportSimulationForML(const SplashScreenResult &params) {
        
    // Ensure output directory exists.
    const std::string outputDir = "output";
    if (!stdfs::exists(outputDir)) {
        if (!stdfs::create_directories(outputDir)) {
            std::cerr << "Error creating output directory: " << outputDir << std::endl;
            return;
        }
    }
    std::string dateTime = getCurrentDateTimeString();
    std::string outputFile = outputDir + "/simulation_export_ml_" + dateTime + ".csv";
    std::ofstream outfile(outputFile);
    if (!outfile) {
        std::cerr << "Error opening output file: " << outputFile << std::endl;
        return;
    }
    
    // Write a CSV header.
    // Header format: run,time,energy,p0_x,p0_y,p0_radius,p1_x,p1_y,p1_radius, ... for all particles.
    outfile << "run,time,energy,bool_dynamic,bool_viscosity," 
        << "gravity_y,overlapParam,interactionParam,viscosity_coeff,sigma,alpha,"
        << "scenarioObj_x,scenarioObj_y,scenarioObj_z,scenarioShapeIndex,"
        << "scenarioObj_min_x, scenarioObj_max_x, scenarioObj_min_y, scenarioObj_max_y";
    for (int i = 0; i < params.numParticles; i++) {
        outfile << ",p" << i << "_x,p" << i << "_y,p" << i << "_radius,p" << i << "_mass";
    }
    outfile << "\n";

    const int   MAX_SEED_ATTEMPTS = 50;          // per run
    const bool  ABORT_IF_STUCK    = true;        // stop whole export after too many failures
    //------------------------------------------------------------------

    for (int run = 0; run < params.numSimulations; /* run++ only on success! */)
    {
        // Create a TemplateApp instance using the splash screen parameters.
        TemplateApp app(params);
        
        // Override simulation parameters with export-mode values.
        app.sim.timeStep = params.timeStep;
        app.sim.dynamic = params.dynamic;
        app.dynamic_convergence_threshold = pow(10.0, params.exponent_convergence_threshold);
        app.sim.viscosity = params.viscosity;
        app.sim.numParticles = params.numParticles;
        app.sim.V0             = params.V0;
        app.sim.L              = params.L;
        app.sim.densify = params.densify;
        app.maxIter = params.maxIter;
        app.sim.use3D = false;  
        app.sim.endTime = params.simulationTime;

        // Set up scenario.
        app.sim.scenarioObjects.clear();
        
        if (params.scenarioShapeIndex == 0) {
            // Circle
            app.sim.scenarioObjects.push_back(
                std::make_unique<Circle>(1.1f, 64, Vector3F(0.0f, 0.0f, 0.0f))
            );
        }
        else if (params.scenarioShapeIndex == 1) {
            // Square
            app.sim.scenarioObjects.push_back(
                std::make_unique<Square>(1.1f, 0.9f, Vector3F(0.0f, 0.0f, 0.0f))
            );
        }
        else /* = 2: Tunnel */ {
            // Tunnel2D: halfLength = 2.0f (visual), halfWidth = sim.L/2
            float halfLen = app.sim.L * 0.5f;
            float halfWid = app.sim.L * 0.5f;
            app.sim.scenarioObjects.push_back(
                std::make_unique<Tunnel2D>(halfLen, halfWid, Vector3F(0.0f, 0.0f, 0.0f))
            );
            app.sim.experiment = Simulation::Experiment::ShearFlow;
            app.sim.periodicX    = params.periodicX;
        }
        //------------------------------------------------------------------
        // 0.1  Try to seed particles --------------------------------------
        //------------------------------------------------------------------
        int tries = 0;
        do {
            app.sim.buildGridDataStructure();
            app.sim.particles2D = app.sim.createRandomParticles2D();
            ++tries;

            if (tries >= MAX_SEED_ATTEMPTS && ABORT_IF_STUCK) {
                std::cerr << "[Export] Could not place particles after "
                        << tries << " attempts – aborting export.\n";
                outfile.close();
                return;
            }
        }
        while (app.sim.particles2D.empty());

        // 0.2  A valid set – continue with normal initialisation ----------
        //------------------------------------------------------------------
        for (auto &p : app.sim.particles2D) p.ix = 0;
        app.reinitializeGlobalState();
        app.sim.minParticles();
        app.sim.buildMassMatrix(app.sim.M);
        app.sim.insertParticlesIntoGrid();
        app.reinitializeGlobalState();

        //------------------------------------------------------------------
        // 1.  Time-stepping loop (exactly what you already have)
        //------------------------------------------------------------------
        const double outputInterval   = 0.01;                 // s
        const int    outputEverySteps = static_cast<int>(
                                        std::round(outputInterval /
                                                    app.sim.timeStep));   // e.g. 1
        int step = 0;
        double t = 0.0;
        app.optimize = true;
        // Simulation loop for one run.
        while (t < app.sim.endTime) {
            //------------------------------------------------------------------
            // 1. Optimization step
            //------------------------------------------------------------------
            Optimization::OptimizationStatus status =
                app.sim.dynamic ? app.energyMinimizationStepDyn()
                                : app.energyMinimizationStep();
        
            // Stop this run if the solver failed to converge.
            if (app.optimize == false) {
                std::cout << "Run " << run
                          << ", time " << t
                          << ": max iterations reached, skipping to next run.\n";
                break;
            }
        
            //------------------------------------------------------------------
            // 2. Advance the simulation state
            //------------------------------------------------------------------
            app.sim.updateScenarioAnimation(app.sim.timeStep);
            t += app.sim.timeStep;                    // <-- increment the clock!
            ++step;
            //------------------------------------------------------------------
            // 3. Diagnostics & CSV output
            //------------------------------------------------------------------
            F energy;
            app.sim.compute_energy(energy);
        
            // Use fmod for a robust “every 0.01 s” check with floating point.
            if (step % outputEverySteps == 0) {
                outfile << run                      << ',' << t            << ','
                        << energy                   << ',' << app.sim.dynamic          << ','
                        << app.sim.viscosity        << ',' << app.sim.gravity(1)       << ','
                        << app.sim.overlapParam     << ',' << app.sim.interactionParam << ','
                        << app.sim.viscosityCoeff   << ',' << app.sim.kernelSigma      << ','
                        << app.sim.alpha            << ',';
        
                // Scenario object pose (assume at least one object exists).
                const auto &obj = *app.sim.scenarioObjects[0];
                outfile << obj.position(0) << ',' << obj.position(1) << ',' << obj.position(2) << ','
                        << params.scenarioShapeIndex << ',';
        
                // Bounding box.
                const BoundingBox bb = obj.getBoundingBox();
                outfile << bb.min_x << ',' << bb.max_x << ',' << bb.min_y << ',' << bb.max_y;
        
                // Particle data.
                for (int i = 0; i < app.sim.numParticles; ++i) {
                    const auto &p = app.sim.particles2D[i];
                    outfile << ',' << p.pos[0] << ',' << p.pos[1] << ',' << p.radius << ',' << p.mass;
                }
                outfile << '\n';
            }
        }
        ++run;
        app.optimize = false;
    }
    
    outfile.close();
    std::cout << "Simulation ML data exported to " << outputFile << std::endl;
}


// void exportSimulationForMLHD5(const SplashScreenResult &params) {
//     std::cout << "Starting exportSimulationForMLHD5" << std::endl;
    
//     // Create the application instance.
//     TemplateApp app(params);
//     app.sim.timeStep = params.timeStep;
//     app.sim.dynamic = params.dynamic;
//     app.dynamic_convergence_threshold = std::pow(10.0, params.exponent_convergence_threshold);
//     app.sim.viscosity = params.viscosity;
//     app.sim.use3D = false;
//     app.sim.endTime = params.simulationTime;
    
//     std::cout << "Simulation parameters set:" << std::endl;
//     std::cout << "  timeStep: " << app.sim.timeStep 
//               << ", dynamic: " << app.sim.dynamic 
//               << ", viscosity: " << app.sim.viscosity 
//               << ", endTime: " << app.sim.endTime << std::endl;
    
//     // Setup scenario.
//     app.sim.scenarioObjects.clear();
//     if (params.scenarioShapeIndex == 0) {
//         std::cout << "Using Circle scenario" << std::endl;
//         app.sim.scenarioObjects.push_back(std::make_unique<Circle>(1.1f, 64, Vector3F(0.0f, 0.0f, 0.0f)));
//     } else {
//         std::cout << "Using Square scenario" << std::endl;
//         app.sim.scenarioObjects.push_back(std::make_unique<Square>(1.1f, 0.9f, Vector3F(0.0f, 0.0f, 0.0f)));
//     }
    
//     // Ensure output directory exists.
//     const std::string outputDir = "output";
//     if (!stdfs::exists(outputDir)) {
//         if (!stdfs::create_directories(outputDir)) {
//             std::cerr << "Error creating output directory: " << outputDir << std::endl;
//             return;
//         }
//         std::cout << "Created output directory: " << outputDir << std::endl;
//     }
//     std::string dateTime = getCurrentDateTimeString();
//     std::string outputFile = outputDir + "/simulation_export_ml_" + dateTime + ".h5";
//     std::cout << "Output file: " << outputFile << std::endl;
    
//     try {
//         // Create the HDF5 file (using Truncate mode)
//         H5::H5File file(outputFile, H5F_ACC_TRUNC);
//         std::cout << "HDF5 file created." << std::endl;
        
//         // ----- Global Dataset Setup -----
//         const int global_cols = 15;
//         const hsize_t maxRows = 10000;
//         // Start with zero rows; extendible along the first dimension.
//         hsize_t global_init_dims[2] = {0, static_cast<hsize_t>(global_cols)};
//         hsize_t global_max_dims[2]  = {maxRows, static_cast<hsize_t>(global_cols)};
//         H5::DataSpace global_dataspace(2, global_init_dims, global_max_dims);
        
//         // Set up chunking.
//         H5::DSetCreatPropList global_prop;
//         hsize_t global_chunk_dims[2] = {1, static_cast<hsize_t>(global_cols)};
//         global_prop.setChunk(2, global_chunk_dims);
        
//         // Create the global dataset.
//         H5::DataSet global_dataset = file.createDataSet("global",
//                                                         H5::PredType::NATIVE_DOUBLE,
//                                                         global_dataspace,
//                                                         global_prop);
//         std::cout << "Global dataset created." << std::endl;
        
//         // ----- Particles Dataset Setup -----
//         int numParticles = app.sim.numParticles;
//         const int particle_cols = numParticles * 4;
//         hsize_t particle_init_dims[2] = {0, static_cast<hsize_t>(particle_cols)};
//         hsize_t particle_max_dims[2]  = {maxRows, static_cast<hsize_t>(particle_cols)};
//         H5::DataSpace particle_dataspace(2, particle_init_dims, particle_max_dims);
        
//         H5::DSetCreatPropList particle_prop;
//         hsize_t particle_chunk_dims[2] = {1, static_cast<hsize_t>(particle_cols)};
//         particle_prop.setChunk(2, particle_chunk_dims);
        
//         H5::DataSet particle_dataset = file.createDataSet("particles",
//                                                           H5::PredType::NATIVE_DOUBLE,
//                                                           particle_dataspace,
//                                                           particle_prop);
//         std::cout << "Particles dataset created." << std::endl;
        
//         // --- Accumulate data from all runs ---
//         std::vector<double> allGlobalData;   // Each row: global_cols elements.
//         std::vector<double> allParticleData;   // Each row: particle_cols elements.
//         int totalRows = 0;
        
//         for (int run = 0; run < params.numSimulations; run++) {
//             std::cout << "Starting simulation run " << run << std::endl;
            
//             // Prepare simulation for this run.
//             app.sim.buildGridDataStructure();
//             app.sim.particles2D = app.sim.createRandomParticles2D();
//             app.sim.insertParticlesIntoGrid();
//             app.reinitializeGlobalState();
            
//             int runRows = 0;
//             double t = 0.0;
//             while (t < app.sim.endTime) {
//                 if (app.sim.dynamic)
//                     app.energyMinimizationStepDyn();
//                 else
//                     app.energyMinimizationStep();
//                 app.sim.updateScenarioAnimation(app.sim.timeStep);
                
//                 double energy;
//                 app.sim.compute_energy(energy);
                
//                 // Prepare global row.
//                 std::vector<double> global_row = {
//                     static_cast<double>(run), t, energy,
//                     app.sim.dynamic ? 1.0 : 0.0,
//                     app.sim.viscosity ? 1.0 : 0.0,
//                     app.sim.gravity(1),
//                     app.sim.overlapParam,
//                     app.sim.interactionParam,
//                     app.sim.viscosity_coeff,
//                     app.sim.sigma,
//                     app.sim.alpha,
//                     app.sim.scenarioObjects[0]->position(0),
//                     app.sim.scenarioObjects[0]->position(1),
//                     app.sim.scenarioObjects[0]->position(2),
//                     static_cast<double>(params.scenarioShapeIndex)
//                 };
//                 allGlobalData.insert(allGlobalData.end(), global_row.begin(), global_row.end());
                
//                 // Prepare particle row.
//                 std::vector<double> particle_row;
//                 for (int i = 0; i < numParticles; i++) {
//                     const auto &p = app.sim.particles2D[i];
//                     particle_row.push_back(p.pos[0]);
//                     particle_row.push_back(p.pos[1]);
//                     particle_row.push_back(p.pos[2]);
//                     particle_row.push_back(p.radius);
//                 }
//                 allParticleData.insert(allParticleData.end(), particle_row.begin(), particle_row.end());
                
//                 runRows++;
//                 t += app.sim.timeStep;
//             }
//             std::cout << "Run " << run << " completed with " << runRows << " rows." << std::endl;
//             totalRows += runRows;
//         }
        
//         // Check accumulated sizes.
//         size_t expectedGlobalSize = totalRows * global_cols;
//         size_t expectedParticleSize = totalRows * particle_cols;
//         std::cout << "Total global data size: " << allGlobalData.size() 
//                   << " (expected " << expectedGlobalSize << ")" << std::endl;
//         std::cout << "Total particle data size: " << allParticleData.size() 
//                   << " (expected " << expectedParticleSize << ")" << std::endl;
//         if (allGlobalData.size() != expectedGlobalSize ||
//             allParticleData.size() != expectedParticleSize) {
//             std::cerr << "Data size mismatch; aborting write." << std::endl;
//             return;
//         }
        
//         // --- Write all global data ---
//         {
//             // Extend the global dataset to accommodate totalRows.
//             hsize_t new_dims[2] = {static_cast<hsize_t>(totalRows), static_cast<hsize_t>(global_cols)};
//             global_dataset.extend(new_dims);
//             H5::DataSpace file_space = global_dataset.getSpace();
            
//             // Select hyperslab starting at (0,0).
//             hsize_t offset[2] = {0, 0};
//             hsize_t count[2]  = {static_cast<hsize_t>(totalRows), static_cast<hsize_t>(global_cols)};
//             file_space.selectHyperslab(H5S_SELECT_SET, count, offset);
            
//             // Create memory dataspace.
//             H5::DataSpace mem_space(2, count);
//             std::cout << "Writing global data (" << totalRows << " rows)..." << std::endl;
//             global_dataset.write(allGlobalData.data(), H5::PredType::NATIVE_DOUBLE, mem_space, file_space);
//             std::cout << "Global data written." << std::endl;
//         }
        
//         // --- Write all particle data ---
//         {
//             hsize_t particle_new_dims[2] = {static_cast<hsize_t>(totalRows), static_cast<hsize_t>(particle_cols)};
//             particle_dataset.extend(particle_new_dims);
//             H5::DataSpace particle_file_space = particle_dataset.getSpace();
            
//             hsize_t offset[2] = {0, 0};
//             hsize_t particle_count[2]  = {static_cast<hsize_t>(totalRows), static_cast<hsize_t>(particle_cols)};
//             particle_file_space.selectHyperslab(H5S_SELECT_SET, particle_count, offset);
            
//             H5::DataSpace particle_mem_space(2, particle_count);
//             std::cout << "Writing particle data (" << totalRows << " rows)..." << std::endl;
//             particle_dataset.write(allParticleData.data(), H5::PredType::NATIVE_DOUBLE, particle_mem_space, particle_file_space);
//             std::cout << "Particle data written." << std::endl;
//         }
        
//         file.close();
//         std::cout << "Simulation ML data exported to " << outputFile << std::endl;
//     } catch (H5::Exception &e) {
//         std::cerr << "HDF5 Exception: " << e.getCDetailMsg() << std::endl;
//     }
// }



