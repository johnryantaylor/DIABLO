# DIABLO
DNS In A Box...Laptop Optimized!

A brief description of DIABLO

diablo.f        Main program, define grid, read/write flow, statistics
input.dat       Set the desired flow parameters
fft.f           Wrapper for calls to FFTW
mpi.f 		 Subroutines for MPI parallelization
grid_def        Define the dimensions of the grid
header          Defines all variables and common blocks

Timestepping subroutines:
periodic.f      Periodic (spectral) in all three directions (done)
channel.f       Periodic (spectral) in x and z, FD in y (done)
duct.f          Periodic (spectral) in x, FD in z and y (not done)
cavity.f        FD in x, y, and z (not done)

In order to run a new case do the following steps:
1.  Create a new run directory
2.  Copy the grid* and input* from an existing case directory and edit them as desired
3.  In the diablo/pre_process directory run create_grid_*.m in matlab or octave to create a new grid for each direction using finite differences (non periodic) and move the created *grid.txt files to the new case directory.
4.  Inside the appropriate timestepping subroutine, edit the create_flow_* subroutine to have the desired initial conditions.
5.  Edit the script “setup_run” and set the rundir variable to match the new directory.
6.  Run “setup_run”.  This will create an executable called “diablo” in your run directory.
7.  Run the code from inside the run directory.  For example, if using MPI, you could do “mpirun -np diablo >output.dat &” which will run a job in the background and write the output to the file output.dat
7.  Some post processing routines can be found in the post_processing directory along with descriptions on their use.

In finite difference (wall bounded) directions, the grid is defined as follows:
(Here the y-direction will be shown, but others should be treated the same)
The points outside the domain are ghost cells and are used to apply the desired boundary conditions.  In the case of one wall-bounded direction, the X and Z velocities (U1 and U3) are defined at the GYF points along with the pressure.  Vertical velocity (U2) is defined at GY points. The GYF points are by definition located halfway between neighboring GY points.

      
 GYF(NY+1)         ---     
                    *            GY(NY+1)
 GYF(NY)      -----Wall-----   
                    *            GY(NY)
 GYF(NY-1)         ---
                    *            GY(NY-1) 
 GYF(NY-2)         --- 
  


 GYF(j+1)          ---
                    *            GY(j+1)
 GYF(j)            ---
                    *            GY(j)
 GYF(j-1)          ---
                    *            GY(j-1)



 GYF(2)            ---
                    *            GY(2)
 GYF(1)      -----Wall-----
                    *            GY(1)
 GYF(0)            ---


Enjoy! 

John

October, 2005  

