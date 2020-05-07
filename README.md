# MATLAB Programming for Finite Element Methods

## Arrangement of the ongoing Toolbox

 We shall establish an iFEM-like package or a simplified version with certain extensions, named mFEM toolbox.

- The toolbox has two important folders：fem and variational.

  - fem: It contains all source files for solving various partial differential equations.
  - variational: A new feature is the variational formulation based programming. It extends the applicaton in fem folder.
                 We have already provided P1, P2 and P3 Lagrange elements for one dimensional and two dimensional problems as well as                    the vectorized problems, e.g. the linear elasticity problem.  

- example: All running examples corresponding to fem and variational are placed in the example folder.

- tool: You can find functions involving visulization, boundary setting, mesh generation, and numerical integration and so on.

- pdedata: It provides information of PDE equations associated with examples in example folder. 

- meshdata: It stories mesh data used in all kinds of examples.

- matlabupdate: Some functions in the updated version of matlab are reconstructed with the same input and output. For example, contains.m ---> mycontains.m.

## Display and marking of polygonal meshes

We present some basic functions to visualize the polygonal meshes, including marking of the nodes, elements and (boundary) edges.

## Auxiliary mesh data and setboundary.m

For the convenience of computation, we introduce some auxiliary mesh data. The idea stems from the treatment of triangulation in iFEM, which is generalized to polygonal meshes with certain modifications.  

## 1-D problems

FEM1D.m and main_FEM1D.m introduce FEM programming of one dimensional problems. The assembly of stiffness matrix and load vector is explained in detail.

## 2-D Poisson equation
The source code of solving the 2-D Poisson equation are presented, see Poisson.m, PoissonP1.m, PoissonP2.m and PoissonP3.m.

## Linear elasticity equations

For linear elasticity problems, we give a unified programming framework. Specifically, 
- The entrices of stiffness matrix are analyzed in the vectorized finite element space;
- The assembly is accomplished in the scalar FE space of each component.

We consider three forms of variational problems. 
  - The first and the second are commonly used in linear elastic problems in the form of strain and/or stress tensors. 
  - The third is just a variant of the second one with Laplacian replacing the strain tensors.
  - For each formulation, sparse assembling indices and detailed explanation are given.

For the first two formulations,  “Neumann”  boundary conditions and Dirichlet boundary conditions are applied. 
For the third fomulation, only Dirichlet conditions are used in view of the practical problems.

## Plate bending problems

- The plate bending problem is to describe the small transverse deformation of thin plates. A special case of the equilibrium equation is the biharmonic equation, which is a typical fourth-order partial differential equation. The nonconforming element is usually applied to reduce the degrees of freedom. Among the nonconforming finite elements in two dimensional case, the Morley element is perhaps the most interesting one. It has the least number of degrees of freedom on each element for fourth order boundary value problems as its shape function space consists of only quadratic polynomials.

    It should be noted that the normal derivative values at the midpoint of interior edge sharing by two triangles have different signs. Apparently, this feature will be inherited by the corresponding local nodal basis functions given by the global ones restricted to the adjacent elements. The problem can be easily resolved by using signed edges.
    
    A more general method is added by introducing signed stiffness matrix and signed load vector. Such prodecure can be extended to other problems with d.o.f.s varing in directions and can be furter applied to polygonal meshes.
	
- Besides Morley element, Zienkiewicz element and Adini element are two other commonly used nonconforming elements. 
  The former is incomplete cubic triangular element and the latter is incomplete bicubic rectangular element.
  In addition to directional problems, all three non-conforming elements (and conforming elements) can be programmed in the unified framework given in the document.

## Mixed FEM
  
   The mixed FEM is applied to solve the biharmonic equation, a special case of plate bending problems.
   

## Variational formulation based programming

  - A variational formulation based programming is shown for 1-D and 2-D problems in Folder variational. The arrangement is entirely       process-oriented and thus is easy to understand. 
  
  - As in FreeFem++, fundamental functions --- int2d.m and int1d.m are designed to match the variational formulation of the underlying    PDEs. These two functions can resolve both scalar and vector equations.
  
  - At present, only Lagrange elements of order up to three are provided, including 1-D problems, Poisson equation, linear elasticity problem, mixed FEM of biharmonic equation and Stokes problem. 
  
  - We remark that the current design can be adapted to find FEM solutions of most of the PDE problems.

  
## Adaptive finite element method and Newest-node bisection

The adpative finite element method (AFEM) is introduced for the Poisson equation with homogeneous Dirichlet 
   conditions.  Each step was explained in detail, viz. the loops of the form: 

           SOLVE -> ESTIMATE -> MARK -> REFINE

The newest-node bisection for the local mesh refinement was clearly stated thanks to the smart idea in iFEM.  
The MATLAB codes are in Folder afem. 

## Multigrid V-cycle method for linear elements

  - We present a multigrid method by adding two grids one by one, so that the pseudo-code of V-cycle method can be obtained directly.
  
  - The MG method can also be analyzed by using subspace correction method proposed by Xu Jinchao (for example, in iFEM). 
    However, it may be more straightforward by adding two grids.

  - The programming of one-dimensional problems is described in detail, and the multigrid program is universal to all linear element problems.
  
  - For 2D and 3D linear elements, only slight changes of 1D problems are needed since they can be regarded as 1D problems.
    See the document for details (two-dimensional problem).


Undo: 

	   
           - mesh generation   
   
           - 3-D FEM      

           - mixed FEM (e.g. Stokes equation)   
  
            - time-dependent problems     

            - ... ...
	    
	    
	    {
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {
    "slideshow": {
     "slide_type": "slide"
    }
   },
   "source": [
    "# Lecture 22: Multi-Layer Neural Network I\n",
    "\n",
    "\n",
    "----\n",
    "\n",
    "### A supervised learning problem\n",
    "We have access to labeled training samples $(\\mathbf{x}^{(i)},y^{(i)})$. Neural networks give a way of defining a complex, highly non-linear form of hypotheses model function $h(\\mathbf{x}; W, b)$, with parameters $W$ and $b$ that we can fit to our data. This nonlinear function is capable of approximating some of the most obscure relations in real life, if we have enough parameters.\n",
    "\n",
    "#### References:\n",
    "\n",
    "* [Simple MNIST numpy from scratch](https://www.kaggle.com/scaomath/simple-mnist-numpy-from-scratch)\n",
    "* [Stanford Deep Learning Tutorial in Matlab](http://ufldl.stanford.edu/tutorial/supervised/MultiLayerNeuralNetworks/)\n",
    "* [3Blue1Brown's video series on Deep Learning](https://www.youtube.com/playlist?list=PLZHQObOWTQDNU6R1_67000Dx_ZCJB-3pi)\n",
    "* [A visual proof that neural nets can compute any function](http://neuralnetworksanddeeplearning.com/chap4.html)\n",
    "* [Backpropagation and chain rule](http://colah.github.io/posts/2015-08-Backprop/)\n",
    "* [Chapter 6: Multiplayer Feedforward Neural Network in *Deep Learning Book* by Goodfellow et al](http://www.deeplearningbook.org/contents/mlp.html)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {
    "slideshow": {
     "slide_type": "slide"
    }
   },
   "source": [
    "## A feedforward neural network\n",
    "\n",
    "Below is one example of a feedforward neural network, the name comes from the fact that the connectivity graph does not have any directed loops or cycles.\n",
    "\n",
    "<img src=\"https://faculty.sites.uci.edu/shuhaocao/files/2019/03/neural_net.png\" alt=\"drawing\" width=\"1000\"/>"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {
    "slideshow": {
     "slide_type": "slide"
    }
   },
   "source": [
    "## How a single neuron works in the $l$-th layer\n",
    "\n",
    "<img src=\"https://faculty.sites.uci.edu/shuhaocao/files/2019/03/neuron-1.png\" alt=\"drawing\" width=\"900\"/>\n",
    "\n",
    "* $\\mathbf{w}$ and $b$: weights and bias\n",
    "* $\\mathbf{a} = (a_1, a_2, a_3)$: input (outputs/activations from the previous layer)\n",
    "\n",
    "This \"neuron\" is a computational unit/node that takes an input $\\mathbf{a}$, and outputs the \n",
    "model function $h^{\\text{single} }(\\mathbf{a}; \\mathbf{w}, b)$ (aka activation):\n",
    "$$ h^{\\text{single} }(\\mathbf{a}; \\mathbf{w}, b) = f(\\mathbf{w}^{\\top} \\mathbf{a} + b) = f\\Big(\\sum_{i=1}^3 w_i a_i +b\\Big) \n",
    "$$\n",
    "The $f(\\cdot)$ is called an \"activation function\", common choices are $\\tanh$, Sigmoid and ReLU:\n",
    "$$\n",
    "\\text{ReLU} (x) = \\max(0, x), \\; \\sigma(x) = \\frac{1}{1 + e^{-x}}, \\; \\tanh(x) = \\frac{e^{x}-e^{-x}}{e^{x}+e^{-x}}\n",
    "$$"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Putting neurons together\n",
    "\n",
    "A neural network is put together by hooking many of our simple \"neurons\", so that the output of a neuron can be the input of another. For example, here is a small neural network (a slice of a bigger network):\n",
    "\n",
    "<img src=\"https://faculty.sites.uci.edu/shuhaocao/files/2019/04/neural_net_3l.png\" alt=\"drawing\" width=\"800\"/>\n",
    "\n",
    "In this figure, we have used circles to also denote the inputs to the network. The circles labeled \"+1\" are called bias units. Layer 1 is called the input layer, and Layer 3 is the output layer (which, in this example, has only one node). The middle layer, Layer 2, is called a hidden layer, because its values are not observed in the training set. \n",
    "\n",
    "The neural network in our example has 2 *input units* (not counting the bias unit), 3 *hidden units*, and *1 output unit*."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Parameters and the forward pass\n",
    "\n",
    "Our neural network has parameters $(W, b) := \\big(W^{(1)},b^{(1)},W^{(2)},b^{(2)}\\big)$.\n",
    "\n",
    "* $W^{(l)} = \\big(w^{(l)}_{ij}\\big)$ to denote the weight matrix, where the entry-$ij$ is associated with the connection between unit $j$ in layer $l$, and unit $i$ in layer $l+1$. Note the order of the indices, $j$ is the closer to the input that this matrix is acting on \n",
    "\n",
    "* $b^{(l)}_i$ is the bias associated with unit $i$ in layer $l+1$. \n",
    "\n",
    "In our example above, we have $W^{(1)}\\in \\mathbb{R}^{3×2}$, and $W^{(2)}\\in \\mathbb{R}^{1×3}$. Note that bias units do not have inputs or connections going into them, we write their output the value `+1` for convenience. When we count the number of units in layer $l$, we do not count the bias unit."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Activation and function compositions\n",
    "\n",
    "We will write $a^{(l)}_i$ to denote the activation or output value of unit $i$ in layer $l$. For $l=1$, $a^{(1)}_i= x_i$ denotes the $i$-th input to this network. Given a fixed set of parameters $(W,b)$, and the input $\\mathbf{x}$, the neural network above defines a model function $h(\\mathbf{x}; W, b)$ made of layers of function compositions that outputs a real number. Specifically, the computation that this neural network represents is given by:\n",
    "\n",
    "$$\\begin{aligned}\n",
    "a_1^{(2)} &= f\\big(w_{11}^{(1)}x_1 + w_{12}^{(1)} x_2  + b_1^{(1)}\\big)  \\\\\n",
    "a_2^{(2)} &= f\\big(w_{21}^{(1)}x_1 + w_{22}^{(1)} x_2 + b_2^{(1)}\\big)  \\\\\n",
    "a_3^{(2)} &= f\\big(w_{31}^{(1)}x_1 + w_{32}^{(1)} x_2 + b_3^{(1)}\\big)  \\\\\n",
    "h(\\mathbf{x}; W,b) &=a =  a_1^{(3)} =  f\\big(w_{11}^{(2)} a_1^{(2)} + w_{12}^{(2)} a_2^{(2)} + w_{13}^{(2)} a_3^{(2)} + b_1^{(2)}\\big) \n",
    "\\end{aligned}\n",
    "$$\n",
    "\n",
    "----\n",
    "\n",
    "## Compact notation: forward pass\n",
    "\n",
    "If we allow the activation function $f(\\cdot)$ to act on vectors in an element-wise fashion: $f([\\mathbf{z}_1,\\mathbf{z}_2,\\mathbf{z}_3])=[f(\\mathbf{z}_1),f(\\mathbf{z}_3),f(\\mathbf{z}_3)]$, then we can write the equations above more compactly as:\n",
    "$$\\begin{aligned}\n",
    "\\mathbf{z}^{(2)} &= W^{(1)} \\mathbf{x} + \\mathbf{b}^{(1)} \\\\\n",
    "\\mathbf{a}^{(2)} &= f(\\mathbf{z}^{(2)}) \\\\\n",
    "\\mathbf{z}^{(3)} &= W^{(2)} \\mathbf{a}^{(2)} + \\mathbf{b}^{(2)} \\\\\n",
    "h(\\mathbf{x}; W, b) &= \\mathbf{a}^{(3)} = f(\\mathbf{z}^{(3)})\n",
    "\\end{aligned}\n",
    "$$\n",
    "More generally, recalling that $\\mathbf{a}^{(1)}=\\mathbf{x}$ also denotes the values from the input layer, then given layer $l$'s activations $\\mathbf{a}^{(l)}$, we can compute layer $(l+1)$'s activations $\\mathbf{a}^{(l+1)}$ as:\n",
    "$$\n",
    "\\begin{aligned}\n",
    "\\mathbf{z}^{(l+1)} &= W^{(l)} \\mathbf{a}^{(l)} + \\mathbf{b}^{(l)}   \\\\\n",
    "\\mathbf{a}^{(l+1)} &= f(\\mathbf{z}^{(l+1)})\n",
    "\\end{aligned}\n",
    "$$\n",
    "By organizing the parameters in matrices and using matrix-vector operations, we can take advantage of fast linear algebra routines to quickly perform calculations in our network."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Loss function for regression\n",
    "\n",
    "Suppose we have a fixed and labeled training set $\\{ (\\mathbf{x}^{(1)}, y^{(1)}), \\dots, (\\mathbf{x}^{(N)}, y^{(N)}) \\}$ of $N$ training examples. For a single training sample and its target value $(\\mathbf{x}, y)$, we define the sample loss function with respect to this single example to be:\n",
    "$$\n",
    "J(W,b; \\mathbf{x},y) = \\frac{1}{2} \\big| h(\\mathbf{x}; W,b) - y \\big|^2,\n",
    "$$\n",
    "or if the label is a vector, \n",
    "$$\n",
    "J(W,b; \\mathbf{x},y) = \\frac{1}{2} \\big\\| h(\\mathbf{x}; W,b) - y \\big\\|^2,\n",
    "$$\n",
    "\n",
    "Then the overall loss function is the mean of the sample losses, plus the regularization term (aka a weight decay term) that tends to decrease the magnitude of the weights $w_{ij}^{(l)}$ but not the biases, and helps prevent overfitting, lastly the $1/2$ factor is added so that upon taking derivate we can get a nice rounded expression without any factors.\n",
    "$$\n",
    "\\begin{aligned}\n",
    "J(W,b)\n",
    "&=  \\frac{1}{N} \\sum_{i=1}^N J(W,b;\\mathbf{x}^{(i)},y^{(i)}) \n",
    "                       + \\frac{\\epsilon}{2} \\sum_{l=1}^{n_l-1} \\; \\sum_{i=1}^{s_l} \\; \\sum_{j=1}^{s_{l+1}} \\left( w^{(l)}_{ji} \\right)^2\n",
    " \\\\\n",
    "&= \\frac{1}{N} \\sum_{i=1}^N \\left( \\frac{1}{2} \\left\\| h(\\mathbf{x}^{(i)}; W,b) - y^{(i)} \\right\\|^2 \\right) \n",
    "+ \\frac{\\epsilon}{2} \\sum_{l=1}^{n_l-1} \\; \\sum_{i=1}^{s_l} \\; \\sum_{j=1}^{s_{l+1}} \\left( w^{(l)}_{ji} \\right)^2,\n",
    "\\end{aligned}\n",
    "$$\n",
    "where $n_l$ denote the number of layers in the network, and $s_l$ denote the number of nodes in layer $l$ (not counting the bias unit)."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## How a neural net works in action\n",
    "\n",
    "We are gonna perform forward passes for a trained neural net with the input layer having 784 input units (28x28 grayscale images), 1 hidden layer with 256 hidden units (neurons), the output layer has 10 units (each represents a class), the activation is ReLU.\n",
    "\n",
    "----\n",
    "\n",
    "#### Implementation Remark:\n",
    "This cost function above is often used both for classification and for regression problems. For classification, we let $y=0$ or $1$ represent the two class labels (recall that the sigmoid activation function outputs values in $[0,1]$; If we were using a $\\tanh$ activation function, we would instead use $-1$ and $+1$ to denote the labels). For regression problems, we first scale our outputs to ensure that they lie in the $[0,1]$ range or $[−1,1]$ range. Most of the time, rescaling inputs is helpful too."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# first the implementation of a vectorized ReLU\n",
    "def relu(x):\n",
    "    return x*(x>0)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# X is the input, of which the second dimension is 784\n",
    "# if X is the trained samples, X.shape = (60000, 784)\n",
    "# if X is a single testing sample, we should make X's shape to be (1, 784)\n",
    "# W is the weight, which is implemented as a list so that\n",
    "# W[0].shape = (784, 256), W[0] maps the input layer to the hidden layer\n",
    "# W[1].shape = (256, 10), W[1] maps the output from the hidden layer to the output layer (10 classes)\n",
    "# b is the bias in each layer, it is also a list so that\n",
    "# b[0].shape = (256,) and b[0] is the bias in the input layer\n",
    "# b[1].shape = (10,) and b[1] is the bias in the hidden layer\n",
    "def h(X,W,b):\n",
    "    # layer 1 = input layer\n",
    "    a1 = X\n",
    "    # layer 1 (input layer) -> layer 2 (hidden layer)\n",
    "    z2 = np.matmul(X, W[0]) + b[0]\n",
    "    # layer 2 activation\n",
    "    a2 = relu(z2)\n",
    "    # layer 2 (hidden layer) -> layer 3 (output layer)\n",
    "    z3 = np.matmul(a2, W[1]) + b[1]\n",
    "    # output layer activation\n",
    "    output = relu(z3)\n",
    "    return output "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# load the trained weights\n",
    "W = np.load('weights.npz')['weights']\n",
    "b = np.load('weights.npz')['bias']"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "X_test = np.load('mnist_test.npz')['X']/255\n",
    "y_test = np.load('mnist_test.npz')['y'].astype(int)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "y_pred = np.argmax(h(X_test, W, b), axis=1)  # pick the biggest activated output unit's index as our prediction"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "print(\"accuracy is:\", np.mean(y_pred == y_test))"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.2"
  },
  "latex_envs": {
   "LaTeX_envs_menu_present": true,
   "autoclose": true,
   "autocomplete": true,
   "bibliofile": "biblio.bib",
   "cite_by": "apalike",
   "current_citInitial": 1,
   "eqLabelWithNumbers": true,
   "eqNumInitial": 1,
   "hotkeys": {
    "equation": "Ctrl-E",
    "itemize": "Ctrl-I"
   },
   "labels_anchors": false,
   "latex_user_defs": false,
   "report_style_numbering": false,
   "user_envs_cfg": false
  },
  "toc": {
   "base_numbering": 1,
   "nav_menu": {},
   "number_sections": true,
   "sideBar": true,
   "skip_h1_title": false,
   "title_cell": "Table of Contents",
   "title_sidebar": "Contents",
   "toc_cell": false,
   "toc_position": {},
   "toc_section_display": true,
   "toc_window_display": false
  },
  "varInspector": {
   "cols": {
    "lenName": 16,
    "lenType": 16,
    "lenVar": 40
   },
   "kernels_config": {
    "python": {
     "delete_cmd_postfix": "",
     "delete_cmd_prefix": "del ",
     "library": "var_list.py",
     "varRefreshCmd": "print(var_dic_list())"
    },
    "r": {
     "delete_cmd_postfix": ") ",
     "delete_cmd_prefix": "rm(",
     "library": "var_list.r",
     "varRefreshCmd": "cat(var_dic_list()) "
    }
   },
   "types_to_exclude": [
    "module",
    "function",
    "builtin_function_or_method",
    "instance",
    "_Feature"
   ],
   "window_display": false
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}

