Copyright © 2024 Danish Center for Earth System Science
Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
and associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, and modify copies of the Software, and 
to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or 
substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

!=============================================================================
Important: fortran libraries netcdf and lapack are needed to run the model
!=============================================================================

!=============================================================================
How to compile and run in linux system (with lapack and netcdf libraries already installed):
!=============================================================================
gfortran -c Parameters.f90 Dimensions.f90 lib_array.f90 Corg_Data.f90 RK4.f90
gfortran -o run_mod Parameters.o Dimensions.o lib_array.o RK4.o -llapack -lnetcdff
./run_mod
!=============================================================================

After run, a netcdf (Out.nc) file is saved containing the main model variables. Some basic model results
can be plotted running the Matlab's file Plot_Model_Results.m located in Plotting folder.




