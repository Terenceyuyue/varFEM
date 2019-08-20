MATLAB Programming for Finite Element Methods

1. Folder FEM1D introduces FEM programming of one dimensional problems. 
   
    The assembly of stiff matrix and load vector is explained in the PDF document 
	"FEM.pdf" in detail (See chapter 1).
	
2. Folder FEM2D includes the source codes of solving the 2-D Poisson equation.
   See chapter 2 for illustration.
   
3. Folder FEM2D-mesh presents some basic functions to show the polygonal meshes, including 
   marking of the nodes, elements and (boundary) edges.
   See chapter 3 for details.
   
4. For the convenience of computation, we introduce some auxiliary mesh data in Folder FEM2D-auxstructure. 
   The idea stems from the treatment of triangulation in iFEM, which is generalized to polygonal meshes with 
certain modifications.  One can refer to Chapter 4 for full explanation.

5. The adpative finite element method (AFEM) is introduced in Chapter 5 for the Poisson equation with homogeneous Dirichlet's 
   conditions.  Each step was explained in detail, viz. the loops of the form: 

           SOLVE -> ESTIMATE -> MARK -> REFINE

The newest-node bisection for the local mesh refinement was clearly stated  thanks to the smart idea in iFEM.  
The MATLAB codes are in Folder AFEM2D. 


Undo: 
      - high-order FEM

      - mesh generation
      
      - Solver: multigrid methods (mg)
      
      - 3-D FEM
      
      - mixed FEM (e.g. Stokes equation)
      
      - time-dependent problems
      
      - ...
