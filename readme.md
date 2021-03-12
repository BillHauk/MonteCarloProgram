# Simulation to Accompany "A Monte Carlo Study of Growth Regressions"
## *Journal of Economic Growth*, 2009

### General Info

This repository consists of five files -- 3 Stata do files and 2 Stata data files.  In order for the simulation to run, all files must be put in the same working directory.  Running the simulation will also create a few extra files in the directory, which can be deleted as needed one the relevant information has been extracted.

Note:  this simulation was first programmed in the mid-2000s in Stata 7.0.  Hence, the first line of code is "Version 7.0".  Some of the commands (especially related to GMM estimators) has been updated in subsequent versions.

### Changing parameters

The "main" program file is MCprogram.do.  All of the parameters that vary in the published paper can be changed in the first 44 lines of code in this file.  The parameters that can be changed are:

#### "True" structural parameters

**alpha** -- the elasticity of output with respect to physical capital in the Solow model used in the data generating process, set by default to 0.27.
**beta** -- the elasticity of output with respect to human capital in the Solow model used in the data generating process, set by default to 0.27.

#### Simulation numbers

**N** -- the number of countries in a single Monte Carlo sample draw, set by default to 69 (the number of countries for which we had a balanced panel data set).
**drawnum** -- the number of Monte Carlo draws in a simulation, set by default to 1000.

#### Measurement error parameters
