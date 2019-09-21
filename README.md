# MATLAB Programming for Finite Element Methods

# Arrangement of the ongoing Toolbox

 We shall estabilish an iFEM-like or a simplified iFEM package, named mFEM toolbox.

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

线弹性问题做了较多修改。我们给出了三种变分形式的程序。
- 第一种和第二种是线弹性问题比较常用的形式，即应变应力张量形式的变分问题。
  - 第一种采用向量有限元空间编程，给出了 sparse 装配指标以及详细的计算说明。
  - 第二种采用标量法编程，即从分量的角度考虑。给出了 sparse 装配指标以及详细的计算说明。
  - 这两种形式的边界条件都包含 “Neumann” 边界和 Dirichlet 边界。

- 第三种是前面变分形式的变形，双线性形式以 Laplace 算子呈现。它的程序类似第二种，但边界条件通常只考虑 Dirichlet 条件。
    
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
