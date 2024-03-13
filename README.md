# Load Balanced Chemistry Model 

This library provides a load-balanced standard and TDAC chemistry model 
for OpenFOAM v2306. This library works with every OpenFOAM solver that uses
the BasicChemistryModel of OpenFOAM. There are currently no known solvers which
are not supported by this library.

## Compilation & Usage

To compile the library source the OpenFOAM environment and execute the 
Allwmake script. 

### Usage

To select the load-balanced chemistry model, select it as the method in the 
chemistryProperties dictionary located in the constant folder. 

```c++
chemistryType
{
    solver            ode;
    method            LoadBalancedChemistryModel;
    // or for TDAC
    method            LoadBalancedTDAC;
}

// Optional Settings
LoadBalancedCoeffs
{
    updateIter  10; // Update the load balance every 10 time steps
                    // default value is every time step
}

// For the load-balanced TDAC model
LoadBalancedTDACCoeffs
{
    updateIter  10; // Update the load balance every 10 time steps
                    // default value is every time step
}
```




