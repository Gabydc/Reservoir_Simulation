# 
# 1D Elliptic FD & FV Solver 
# Incompressible Flow Solver - Gaby 24/04/17

# Initialization Parameters, Grid, ...
# Grid Configurations
# |---X---|---X---|...
# Note: # of interfaces = # cells +1

L = 1.0;    # Length of the Reservoir [m]N = 20;    
N = 20;     # Number of Grid Cells
PL = 1;     #INPUT BC
PR = 0;     #INPUT BC

DX = L/N;   # Grid size
x =     range(DX/2, L - DX/2, DX);   #Location of Grid centers
xi =    range(0, L, DX); #Location of interfaces