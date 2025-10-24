Notes for running and evaluating the ocean running on a gpu.

Right now we are just doing this for the split-explicit version of the code - in user_nl_mpaso set 
  config_time_integrator = 'split_explicit'

To get a control for testing this:
  git clone -b develop https://github.com/EarthWorksOrg/EarthWorks.git EarthWorks-cont

To get the gpu branch:
  git clone -b develop https://github.com/EarthWorksOrg/EarthWorks.git EarthWorks-gpu
  cd EarthWorks-gpu
  edit .gitmodules and replace the mpas-ocean and mpas-framework entries with

[submodule "mpas-ocean"]
        path = components/mpas-ocean
        #url = https://github.com/EarthWorksOrg/mpas-ocean.git
        url = https://github.com/Pranay-Reddy-Kommera/mpas-ocean.git
        fxDONOTUSEurl = https://github.com/EarthWorksOrg/mpas-ocean.git
        fxrequired = ToplevelRequired
        #fxtag = mpaso-ew2.5.006
        fxtag = mpaso-openacc-integration
        
[submodule "mpas-framework"]
        path = components/mpas-framework
        #url = https://github.com/EarthWorksOrg/mpas-framework.git
        url = https://github.com/Pranay-Reddy-Kommera/mpas-framework.git
        fxDONOTUSEurl = https://github.com/EarthWorksOrg/mpas-framework.git
        fxrequired = ToplevelRequired
        #fxtag = mpasfrwk-ew2.5.000
        fxtag = mpaso-openacc-integration

Evaluating the gpu solution.

We build three runs to evaluate the gpu code. We use the nvhpc compiler.

1) Build and run the control code (this is on a cpu).
   ~dazlich/cpu_o120.csh will build and submit a job

2) Build and run the gpu branch code on cpu.
   modify ~dazlich/cpu_o120.csh - change the runname and change EarthWorks-cont to EarthWorks-gpu, then run the script.

3) Build and run the gpu branch code on gpu.
   ~dazlich/gpu_o120.csh will build and submit a job

Notes on gpu compilation - check the bld/ocn.bldlog file to make sure the ocean code is being compiled with the gpu flags you think you are using.

GPU compilation flags can be modified in mpas-ocean/cime-config/buildlib

Comparing solutions

These runs are all one month simulations at 120km on one derecho node.

First, we want to make sure run1 and run2 are producing the same solution. 

ncdump -v kineticEnergyCellAvg run_1_globalStats_file.nc | tail -12
ncdump -v kineticEnergyCellAvg run_2_globalStats_file.nc | tail -12

These are writing daily and domain averaged values and should agree for most of the decimal places for the entire simulation.

If the gpu run is executing
ncdump -v kineticEnergyCellAvg run_3_globalStats_file.nc | tail -12
This will differ a bit more but we can expect at least 4 or 5 digit agreement for the values at the end of the simulation.

Questions? - donald.dazlich@colostate.edu
