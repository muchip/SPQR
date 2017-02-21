SParse Quadrature Routines (SPQR)

A Package which implements the anisotropic sparse grid quadrature for arbitrary 
downward closed index sets. To provide maximum flexibility, the underlying 
univariate quadrature rules as well as the criterion for the sparse index set can 
be defined by the user. Optimised versions of the total degree index set and the 
hyperbolic cross index set are provided by default. 
Currently, the package comes with a simple Matlab interface.
The theoretical foundation of the anisotropic sparse grid quadrature
can be found in preprint.pdf.

INSTALLATION

For the installation, you will need make, a C++ compiler, Eigen and the Matlab-mex compiler.
Particularly, you will have to set the correct path to your C++ compiler, to your Eigen library 
and to your Matlab-mex compiler. Afterwards, just execute make.
