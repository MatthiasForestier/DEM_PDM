#ifndef SPLASHSCREEN_H
#define SPLASHSCREEN_H

// Displays a splash screen allowing the user to choose between 2D and 3D modes.
// Returns true if 3D mode is selected, false if 2D mode is selected.
struct SplashScreenResult {
    bool use3D;       // true for 3D mode, false for 2D mode
    // int numParticles; // number of particles for the simulation
};

SplashScreenResult showSplashScreen();

#endif // SPLASHSCREEN_H
