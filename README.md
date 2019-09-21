# MATLAB Programming for Finite Element Methods

# Arrangement of the ongoing Toolbox

 We shall establish an iFEM-like or a simplified iFEM package, named mFEM toolbox.

- The toolbox has two important folders：fem and variational.

  - fem: It includes all kinds of source code.
  - variational: A new feature is the variational formulation based programming. It extends the applicaton in fem folder.
                 At present, only one-dimensional problems are provided.  

- example: All running programs are placed in the example folder.

- tool: You can find functions involving visulization, boundary setting, mesh generation, and numerical integration and so on.

- pdedata: It provides information of equations associated with examples in example folder. 

- meshdata: It stories mesh data used in all kinds of examples.

## 1-D problems

Folder FEM1D introduces FEM programming of one dimensional problems.    
    The assembly of stiff matrix and load vector is explained in the PDF document 
	"FEM.pdf" in detail (See chapter 1).

## 2-D Poisson equation
Folder FEM2D includes the source codes of solving the 2-D Poisson equation.
   See chapter 2 for illustration.

## Display and marking of polygonal meshes
Folder FEM2D-mesh presents some basic functions to show the polygonal meshes, including 
   marking of the nodes, elements and (boundary) edges.
   See chapter 3 for details.

## Auxiliary mesh data and setboundary.m
For the convenience of computation, we introduce some auxiliary mesh data in Folder FEM2D-auxstructure. 
   The idea stems from the treatment of triangulation in iFEM, which is generalized to polygonal meshes with 
certain modifications.  One can refer to Chapter 4 for full explanation.

## Adaptive finite element method and Newest-node bisection
The adpative finite element method (AFEM) is introduced in Chapter 5 for the Poisson equation with homogeneous Dirichlet's 
   conditions.  Each step was explained in detail, viz. the loops of the form: 

           SOLVE -> ESTIMATE -> MARK -> REFINE

The newest-node bisection for the local mesh refinement was clearly stated  thanks to the smart idea in iFEM.  
The MATLAB codes are in Folder AFEM2D. 

## Linear elasticity equations

For linear elasticity problems, we give a unified programming framework. Specifically, 
- The entrices of stiffness matrix are analyzed in the vectorized finite element space;
- The assembly is accomplished in the scalar FE space of each component.

We consider three kinds of variational formulation. 
  - The first and the second are commonly used in linear elastic problems in the form of strain and/or stress tensors. 
  - The third is just a variant of the second one with Laplacian replacing the strain tensors.
  - For each formulation, sparse assembling indices and detailed explanation are given.

For the first two formulations,  “Neumann”  boundary conditions and Dirichlet boundary conditions are applied. 
For the third fomulation, only Dirichlet conditions are used in view of the practical problems.
    
## Variational formulation based programming for 1-D problems
We present a variational formulation based programming for 1-D problems in Folder variational1D. The arrangement is entirely process-oriented and thus easy to understand (see Section 1.5). 


PS:

     作者最近忙于弹性力学的总结与程序实现，对前面各章节的重新整理将会在未来的适当时间完成。

Undo: 

           --- Variational formulation based toolbox
   
           - high-order FEM
	   
           - mesh generation   
   
           - Solver: multigrid methods (mg)   
   
           - 3-D FEM      

           - mixed FEM (e.g. Stokes equation)   
  
            - time-dependent problems     

            - ... ...
