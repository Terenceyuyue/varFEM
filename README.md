# MATLAB Programming for Finite Element Methods

# Toolbox 结构安排
之前上传的文档已经删除，本文将重新整理文档，以形成类似 iFEM 的工具箱。

该工具箱有两个重要的文件夹：fem 和 variational。
后者存放 “基于变分形式编程的程序”，目前只给出一维问题的处理。
fem 文件夹存放常规编程的各种函数文件。作者最近在处理弹性力学问题，所以先上传该部分内容。

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
