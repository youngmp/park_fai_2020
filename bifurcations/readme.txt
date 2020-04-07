This folder contains all ODE files and associated AUTO files. 


========================
for these bifurcation diagrams, it is best to use the following numerics options.

Ntst	= 15 (default)
Nmax	= 200 (default)
Npr	= 50 (default)
Ds	= 1e-10
Dsmin	= 1e-10
Ncol	= 4 (default)
Dsmax	= 0.1

adjust Ntst, Nmax, Npr, Par Min, and Par Max as needed. This Dsmax value will make AUTO run multiple passes over the branches, so it may help to set it to 0.001 to prevent multiple passes, but also to increment at a reasonable pace.


allinfo_2par.dat: all plotting information for 2 parameter diagram
*.auto: auto diagram files for XPP
p1=*.dat: bifurcation diagram (allinfo) with fixed p1=*
ze=*.dat: bifurcation diagram (allinfo) with fixed ze=*

p1 stands for phi1, not to be confused with the p1 parameter for force-position in the agent-based simulations

pi4=*_pi5=*.dat: allinfo two paramter bifurcations diagrams with pi4 and pi5 set to the wildcard values. Used to generate fig. 4 panels B, C, D, E. Panel A was generated using code from cusp.py, which detects changes in the number of equilibria using a root finding method (brentq).




========================
force.ode: original ODE file used to generate 2 parameter bifurcation diagrams
The following were not used
#force_saddle.ode: Eliminate 1 parameter assuming a saddle-node bifurcation (failed to work to track cusps).
#force_cusp.ode: Eliminate 2 parameters assuming a cusp (failed to work -- had to turn to root-finding to "track" cusp).

*.auto files are not guaranteed to work.

__init__.py: used for importing this folder as a python module.

========================
