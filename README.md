# DifEqMaster
A personal C++-based framework for solving some time-domain partial differential equation systems 

## Introduction
This project began in 2002 as an attempt to summarize and systematize some multidimensional differential equation solvers that I've made for numerical simulation of vacuum electronic radiation sources. I wanted to put my standalone code for simulating excitation processes in free electron lasers and backward wave oscillators under a universal platform. That C++ Builder based platform helped me with setting dozens of parameters, following the simulation process in realtime and presenting the result graphs in a more unified and convenient way. Differential equation systems for some other types of lasers were added for simulation in the following years with general focus being shifted into solid state laser structures with 1D and 2D distributed feedback. Various tools and techniques were utilized for smoothing both computation process and interpretation of the results. Some notable additions include multithreading, FFT-based spectrum analysis, particle swarm techniques for automatic optimization of laser parameters. Combined results of this project became the basis of my Master and Ph.D. degrees as well as over 20 published articles in peer-reviewed journals and conference speeches from 2004 to 2023.

## Project Structure
### Common Files
The main idea of the project was to share input/output procedures between various differential equation solvers. So, naturally, we have several common parts for input/output:

* MainForm (.h/.cpp)  for putting up the main window and setting up everything
* EquGraph (.h/.cpp) for scaling and drawing graphical results
* EquBaseDump, DumpStorage (.h/.cpp)  for saving dumps

We also have some common files with abstract classes:
* EquBaseMaster, EquBaseThread and EquBaseTask (with .h/.cpp)
  
Task-like classes are where all the real computation takes place, so they were intended to be fast and totally free of input/output. Some of them acquired multi-thread versions when multi-thread capable CPUs became available to me. Every equation system here is expected to be supplied with a Master class for handling its dozens of parameters, setting everything up and doing other once-in-a-while and non-time-critical chores.

Additionally, there are some libraries with mathematical tools:

* FFTUnit (.h/.cpp) with Fourier transformations
* SplineUnit (.h/.cpp)   with spline approximations
* MultiMaxPSO (.h/.cpp)  and other PSO-named files with my implementations of the Particle Swarm Optimization technique
* CephesLib (.h/.cpp)  with Fresnel integral computation functions that I’ve took from Cephes Mathematical Library. It is quite an old library and https://github.com/jeremybarnes/cephes is the only working reference for it that I’ve found
  
### Free electron models

Direct simulation of free electron systems can be complicated by a huge mismatch of scale. Electron beams can propagate for more then 20 meters in some free electron lasers, while it can be critically important to keep track of electron’s positions with sub-nanometer precision. In those systems excitation process can be simulated with a technique known as “particle-in-cell” (PIC) that groups a large number of free electrons in a “super particle” with a macroscopic size, but similar electron phases. This technique separates “approximate positions” of the electron beam parts and their phases, which effectively are “precise positions relative to the electromagnetic wave”.

In this part of the project we implement several systems of 2D and 3D differential equations with partial derivatives to simulate and optimize various free electron devices like free electron lasers and backward wave oscillators. They take into account relativistic effects, Doppler effects, transverse electron oscillations, thermal velocity distribution inside electron beams, channeling of electromagnetic radiation by electron beams. 

All of this is split into following files:

Equ2WaveElModel, Equ2WaveTask, Equ2WaveElModel, Equ2WaveTask, Equ25DEl01, Equ25DEl02, Equ25DEl03, Equ25DEl04, Equ25DEl05, Equ25DEl06, EquOldLBVUnit (all names are pairs of .h/.cpp files) 

### Solid-state distributed feedback laser models

This part of the project describes various laser structures with 1D and 2D distributed feedback. The main goal here is demonstrating the excitation process that begins with spontaneous radiation and results in steady-state radiation regimes.  Simulated and demonstrated effects here include spatial synchronization of radiation across two dimensions, interference of two propagating waves that results in spatial hole burning, synchronization of multichannel laser structures, mutual scattering of several electromagnetic waves.
The complicated part here is producing a mathematically stable procedure for calculating mutually scattering waves. Most typical procedures like Runge-Kutta methods for differential equations proved to be unstable for considered cases. Consequently, we use a semi-implicit procedure for making a time domain step. Additionally, electromagnetic waves that propagate in different directions complicate splitting the step procedure into separate procedures for separate spatial regions for multithreading.

This part of the project is split into following files:

EquOptic2DMasterForm, EquOpticTask, EquOptic2D, EquOptic2D_1DwB, EquOptic2D_1DwB_Mod1, EquOptic2D_MultiMedia, EquOptic2D_X_01, EquOptic2D_X_02, EquOptic2D_X_03, EquOptic_X_04, EquOptic2D_X_MultiMedia, EquOptic2D_X_Multimedia2, EquOptic2D01, EquOptic2D02, EquOptic2D03, EquOptic2D04, EquOpticTask_2BraggsSolidMedia, EquOpticTask_2D1D2DSolidMedia, EquOpticTask_2Way, EquOpticTask_BMB, EquOpticTask_BMB_holes, EquOpticTask_MultiMedia, EquOpticTask_Wide1DBragg, EquPSO2DBraggSolver (all names are pairs of .h/.cpp files) 

## Addendum/Disclamer/Conclusion

As mentioned before, this project became the basis for my Master and Ph.D. degrees works as well as a variety of published scientific articles. Most of the articles can be found here:
https://www.researchgate.net/scientific-contributions/V-R-Baryshev-32973328

With that being said, I must admit that both the project structure and its C++ Builder platform are quite outdated. Partially, this is a result of the project's 20 years lifespan. Something was being altered and added during this whole time and that creates a mess. Also, with me being the sole owner, producer, designer, QA and end user, it was tempting to skip some household cleaning with "I only need to run it once and get the result" attitude. Regrettably, this gave way for some more mess.




