# MATLAB Programming for Finite Element Methods

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
We consider the programming of linear elasticity equations in Chapter 6, which is a classical vector equation arising from 
   elastic mechanics. The assembly of stiff matrix and load vector is clearly presented following the idea of scalar equations. 
   The presented code may be more accepatable than the one in iFEM. We remark that the idea can be directly applied to 
    Stokes equations which will be introduced later for the mixed FEM.  The MATLAB codes are included in Folder elasticity.
    
## Plate bending problems (Morley element)

薄板弯曲问题导出双调和方程，这是典型的四阶问题，通常用非协调元计算，以减少自由度。最简单的非协调板元是 Morley 元，它的自由度最少且刚度矩阵很容易计算。 我们给出了板弯问题的 Morley 元程序。

需要注意的是，Morley 元边中点的法向自由度对相邻两个单元来说是异号的，因而整体节点基限制在相邻单元上会异号。通过将边符号化，这个问题很容易得到解决。

文档中详细介绍了 Morley 元的计算及装配问题 (章节名: Kirchhoff 板弯问题)。

    
## Variational formulation based programming for 1-D problems
We present a variational formulation based programming for 1-D problems in Folder variational1D. The arrangement is entirely process-oriented and thus easy to understand (see Section 1.5). 



Undo: 

           --- Variational formulation based toolbox
   
           - high-order FEM
	   
           - mesh generation   
   
           - Solver: multigrid methods (mg)   
   
           - 3-D FEM      

           - mixed FEM (e.g. Stokes equation)   
  
            - time-dependent problems     

            - ... ...
