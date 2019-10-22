Portable Model Readme

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This is a portable radiative-convective model, which can be run in 1D or 2D (latitude and pressure).

For the first time user, the important file is 'Call Fortran from Python [Interpolation Separate].py'. In this Python script, you can set the values of a large number of model parameters, including but not limited to CO2, H2O, and O3 amounts, cloud altitude, fraction, and optical thickness, lapse rate (NOTE: I've disabled convection in this version to allow you to try the alternative method), surface albedo, and surface pressure.

This Python script uses these parameters to write input files to the local directory, which can then be read by 'wrapper v1.1 2018-08-13 return to backup after VARIABLES module issue.f95'. This is the main Fortran program. It provides a wrapper in which layer and level p, T, z, etc. are calculated and the necessary information is prepared in the correct form to be used by 'rrtm.f', which performs the radiative transfer calculations.

At the end of the run, the model writes an output file stamped with the date and time into a folder called '_Raw Output Data'. The output file contains the final temepratures, fluxes, heating rates, and several other variables.

While the model runs, you'll see printed the tropopause layer and the maximum heating rate above the tropopause. This should approach 0 as equilibrium is reached. Currently it's a bit buggy and often mis-identifies the maximum heating rate, but it still provides an idea of how the run is progressing. Don't worry about this too much - you can see the actual heating rates when you plot the output file.

To plot and/or print data from the output file, drag and drop it into the folder '_Current Output'. Then run the Python script 'Read and plot 2D RCM output data (compare multiple files).py'. This will produce basic plots of T, heating rate, fluxes, and greenhouse gase mixing ratios. If you want to print values of any variables they are available for each file, layer, and latitude column as described in the comments of that file.

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

A quick guide to some of the other files and folders, which you shouldn't need to use but might be curious about:

	- 2D RCM GitHub.dSYM
		- Keeps track of debugging symbols
	- _Useful Data
		- I move output files that I know are useful results here to keep track of them
		- To plot from here, just change the 'directories' list in 'Read and plot 2D RCM output data (compare multiple files).py'
	- Debug 
		- Folder created when compiling in Debug mode; IGNORE
	- Input Data for RCM Interpolation
		- This contains numpy arrays that were created earlier when I interpolated the ERA-Interim data onto a smaller grid that's large enough for the purpose here
		- At the start of the running the model, these arrays are used to create input files of the size needed for the current model parameters
	- Input Distributions
		- This folder contains the files that the main Fortran program reads directly as inputs for clouds, greenhouse gases, and (optionally) temperatures
	- Latitudinal Distributons
		- Latitudinal distributions scanned from figures in other papers
		- Currently only Mason's lapse rates are used from here
	- Release
		- Folder created when compiling; IGNORE
	- src
		- The original Fortran source files needed to compile and run RRTM
	- Earth RCM Parameters
		- This is where 'Call Fortran from Python [Interpolation Separate].py' writes the model parameters
		- This is read directly by the main fortran program
	- MYSUBS.f95
		- A Fortran module where I keep subroutines that I've written to keep the main Fortran program more readable and easier to edit
		- Includes subroutines to print the current time, perforam the convective adjustment, perform SW radiative transfer, create cloud input and output files, and write the main output file
	- VARIABLES.f95
		- A module that contains all the variable declarations
	- 'wrapper v1.1 2018-08-13 return to backup after VARIABLES module issue.f95'
		- The main Fortran program that reads in the input files, calls RRTM to do the radiative transfer, calls the subroutines to get to radiative equilibrium and write output files, and prints progress to the standard output window