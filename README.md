# varFEM: Variational formulation based programming for Finite Element Methods in MATLAB

If you have any problems, please contact me via the email: terenceyuyue@sjtu.edu.cn

[TOC]
 
 --------------

## Intentions

We intend to develop the "variational formulation based programming"  in a similar way of FreeFEM, a high level multiphysics finite element software. The similarity here only refers to the programming style of the main program, not to the internal architecture of the software.


Example: The stiffness matrix for the bilinear form  $\int_{T_h} \nabla u \cdot \nabla v {\rm d}\sigma$ can be computed as follows.

```
  Vh = 'P1';  quadOrder = 5;
  Coef  = 1;
  Test  = 'v.grad';
  Trial = 'u.grad';
  kk = assem2d(Th,Coef,Test,Trial,Vh,quadOrder);
```

## Variational formulation based programming

  - A variational formulation based programming is shown for 1-D, 2-D and 3-D problems. The arrangement is entirely process-oriented and thus is easy to understand. 
  
  - **As in FreeFem, fundamental functions --- int1d.m, int2d.m and int3d.m are designed to match the variational formulation of the underlying PDEs. These functions can resolve both scalar and vectorial equations (see also assem1d.m, assem2d.m and assem3d.m).**
  
  - Only Lagrange elements are provided in the current version. 
  
  - We remark that the current design can be adapted to find FEM solutions of most of the PDE problems.
  
## Examples

  - Poisson equation
  - Linear elasticity problem
  - Biharmonic equation (mixed element)
  - Stokes problem
  - Heat equation
  - Navier-Stokes equation
  - Variational inequalities
  - Eigenvalue problems
  - Unified implementation of adaptive finte element methods: Poisson equation
  - Unified implementation of FEMs involving jumps and averages in bilinear forms: C0 interior penalty method for the biharmonic equation
