# TDGL.jl
Time-dependent Ginzburg-Landau simulations with the use of Julia 1.0.

This software is released under the MIT License, see LICENSE.

```
add https://github.com/cometscome/TDGL.jl
```

# test
You need the output directory. 
Before doing the simulation, 

```
mkdir output
```
is needed.

Then, 
```
using TDGL
TDGL.test2()
```
You can simulate $25 \xi_0 \times 25 \xi_0$ superconducting nano island with the magnetic field $H = 0.1$ at $T = 0.2T_c$.
