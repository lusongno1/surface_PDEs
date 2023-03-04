# surface_PDEs Project
Welcome to My Surface_PDEs Project! This project is designed for computation of partial differential equations and eigenvalue problems on surfaces. 
For an introduction to the calculation method, please refer to [Lu S, Xu X. A geometrically consistent trace finite element method for the Laplace-Beltrami eigenvalue problem [J]. ](https://arxiv.org/abs/2108.02434)

## Installation&Dependencies
To fully run this project, some packages need to be installed, listed below. Please follow the corresponding installation instructions to install these softwares and their respective dependencies.
### [DROPS](https://www.igpm.rwth-aachen.de/forschung/drops)
The key module in DROPS_highquad is forked from [https://github.com/56th/drops](https://github.com/56th/drops), we need to follow the instructions in it to install it.

To build DROPS_highquad, follow these steps:
1. Clone this repository
2. Run `sh cbmaker.sh` to install it
3. Start the project by Code::Blocks
4. Run cpp files

### [PHG](http://lsec.cc.ac.cn/phg/)
The PHG version used is [http://lsec.cc.ac.cn/phg/download/phg-0.9.5-20200521.tar.bz2](http://lsec.cc.ac.cn/phg/download/phg-0.9.5-20200521.tar.bz2), and the high-order numerical quadrature interface of subsequent versions has been fine-tuned. It is a key dependency of DROPS_highquad.


### [NGSXFEM](https://github.com/ngsxfem/ngsxfem)
Toolbox [ngsxfem](ttps://github.com/ngsxfem/ngsxfem) (the packages builds on NGSolve) or [tutorials]( https://github.com/ngsxfem/ngsxfem-jupyter) are needed, you can easily try it out. If you don't want to install it locally you can check it out in the cloud through binder, here: [https://mybinder.org/v2/gh/ngsxfem/ngsxfem-jupyter/HEAD?filepath=tutorials.ipynb](https://mybinder.org/v2/gh/ngsxfem/ngsxfem-jupyter/HEAD?filepath=tutorials.ipynb) , see Section 5 on Surface PDEs in it.

### [MultiParEig](www.mathworks.com/matlabcentral/fileexchange/47844-multipareig)
[MultiParEig](www.mathworks.com/matlabcentral/fileexchange/47844-multipareig) is a software package for computing multiple eigenvalues and corresponding eigenvectors of parameter-dependent matrices. It is developed by researchers from the University of Antwerp and KU Leuven in Belgium.

## Usage
For DROPS_highquad, the main programs we use are all in the `src/surfactant` folder, and the parameter input files used are in the `/param/surfactant/surfactant` folder.
- surfactant.cpp: Computes the solution to the Laplace-Beltrami equation on a surface.
- eigenvalue_problems.cpp: Computes generalized algebraic matrices for the Laplace-Beltrami eigenvalue problems on surfaces.
- eigfun.cpp: Plot the eigenfunction of some operator on the surface.

The MATLAB scripts is mainly used to calculate the generalized algebraic eigenvalue problem by some different methods. The program entry refers to the m-files named beginning with 'test' under the `MATLAB` folder. The subfolders under the `MATLAB` folder store some data of the stiffness matrix, mass matrix calculated by the numerical method from DROPS and NGSXFEM. 

The Python file eig.py needs to be run under the environment of NGSXFEM, which gives the relevant algebraic matrix of the isoparametric finite element method.

For more details about program operation, please contact [lusong@lsec.cc.ac.cn](lusong@lsec.cc.ac.cn) or [xmxu@lsec.cc.ac.cn](lusong@lsec.cc.ac.cn).



## Contributing

If you'd like to contribute to the Project, please follow these guidelines:

1. Fork this repository
2. Create a new branch for your changes
3. Make your changes and commit them
4. Submit a pull request

## License

My Project is licensed under the MIT License. See `LICENSE` for more information.
