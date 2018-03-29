# FES_codes
Codes for finding, manipulating, and pulling structures from LE4PD or any other free-energy surfaces (FESs).

What does what:

print_fes.py -- Simply a python script that takes a processed trajectory (anly_\*.dat) in two dimensions, computes the 2D-histogram of the data, then calculates the FES of the distribution, and plots it as a filled-contour plot and writes it to file.

coords.py -- Plots the FES and then uses the matplotlib.pyplot.ginput function to interactively select and write the selected coordinates to file.

FTS-Ensemble-Evolve.py -- A non-modularized version of the finite-temperature-string (FTS) method for finding the most probable energy pathway along a 2D energy surface. Modeled after the information presented in Weinan, Ren, and Vanden-Eijnden, *J. Phys. Chem. B* **2005**, *109*, 6688-6693.

average_structures.f95 -- FORTRAN code for finding the average conformation of a given ensemble of structures. Probably better to use the GROMACS option from gmx rmsf : gmx rmsf -f structs.pdb -s top.tpr -ox average.pdb -fit no < input.inp

ensemble_finder.f95 -- Deprecated code for finding the ensemble of structures correpsonding to the given coordinates. Likely better and more efficient to use ...

test_ensemble_finder.f95 -- Should probably be re-named, but this version reads the given coordinates and the processed trajectory simultaneously to find the frames that acheive the given coordinates, to within a threshold value. Writes the frames corresponding to each structure to a different .ndx for extraction by pull_structures_by_time.sh

pull_structures_by_time.sh -- Shell script for pulling selected structures from a compressed XTC trajectory, given the desired frames and the number of desired configurations. Calls test_ensemble_finder.f95 to find the ensemble of configurations, then pulls those configurations from the trajectory using gmx trjconv and the respective .ndx files. The average confguration of each ensemble is found using gmx rmsf: gmx rmsf -f structs.pdb -s top.tpr -ox average.pdb -fit no < input.inp. Using the -fit no option is quite important.
**Note**: assumes GROMACS version 5 is being used. Could require tweaking the "mult" variable for other versions of GROMACS.

structure_finder.f95 -- Uses a distance criterion to find the closest structure from a set of structures to the given structure.

interpol_traj.f95 -- Determines the linear interpolation between two given PDB files and writes the interpolated trajectory to a single PDB file.

overlay.py -- Overlays the interpolated trajectory and the FES calculated from print_fes.py. Useful for examining how closely the interpolated trajectory follows the most probable pathway.

Any code present in the directory that has not been described either is either deemed unimportant (to me) or is deprecated.
