# DifEqMaster
A personal C++-based framework for solving some time-domain equation systems

## Introduction
This project began in 2002 as an attempt to summarize and systematize some multidimensional differential equation solvers for that I've made for computing a number of problems in free electronic devices. I wanted to put some standalone code for simulating excitation processes in free electron lasers and backward wave oscillators under a universal platform. That C++ Builder based platform helped me with setting dozens of parameters, following the simulation process in real-time and presenting the result graphs in a more unified and convenient way. Equations for some other laser systems were added for simulation in the following years with general focus being shifted into solid state laser structures with 1D and 2D distributed feedback. Various tools and techniques were utilized for smoothing both computation process and interpretation of the results. Some notable additions include multithreading, FFT-based spectrum analysis, particle swarm techniques for automatic optimization of laser parameters. Combined results of this project became the basis of my 2004 Master and 2010 Ph.D. degrees as well as over 20 published articles in peer-reviewed journals and conference speeches from 2004 to 2023.

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
* MultiMaxPSO (.h/.cpp)  and other PSO-named files with my implementations of the Partical Swarm Optimization technique
* CephesLib (.h/.cpp)  with Fresnel integral computation functions that I’ve took from Cephes Mathematical Library. It is quite an old library and https://github.com/jeremybarnes/cephes is the only working reference on it that I’ve found
  
### Free electron models
Direct simulation of free electron systems can be complicated by a huge mismatch of scale. Electron beams can propagate for more then 20 meters in a free electron laser, while it can be critically important to keep track of electron’s positions with sub-nanometer precision.
Excitation process in free electron lasers, backward wave oscillators and a number of other free electron devices can be simulated with a partical-in-cell (PIC) technique that groups a large amount of free electrons in a “super particle” with a macroscopic size, but similar phase. This technique effectively separates “positions” of the electron beam parts and their phases
described as three wave interaction process between two electromagnetic waves and one electron beam current wave

Classes for free-electron models
