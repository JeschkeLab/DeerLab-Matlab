
# DeerLab has been migrated to Python
You can find the latest version of the source code in the [DeerLab repository](https://github.com/JeschkeLab/DeerLab), where it will be continued to be developed. All DeerLab pre-releases written in MATLAB are still available in this repository.

### Requirements
DeerLab requires the following products:

  * MATLAB (R2017a or newer) (see <https://www.mathworks.com/products/matlab.html>)
 
DeerLab will use the following product if installed:
 
  * Optimization Toolbox (see <https://www.mathworks.com/products/optimization.html>)

### Setup

In order for MATLAB to access the DeerLab functions, the path to the DeerLab installation folder must be set in MATLAB.

**Option 1:** Add DeerLab path via MATLAB's IDE

1) On the ``Home`` tab, in the ``Environment`` section, click ``Set Path``. 

2) Click ``Add with Subfolders...`` and select the ``DeerLab\functions`` directory. 

3) Click ``Save`` to save the current MATLAB search path and exit via ``Close``.

**Option2:**  Add DeerLab path at startup

1) Open (or create) the ``startup.m`` file in the default ``\MATLAB`` directory.

2) Add the following lines of code:

       addpath('mypath/DeerLab/functions')

3) Save ``startup.m`` and restart MATLAB.

### Citation

A publication about DeerLab is [available here](https://doi.org/10.5194/mr-2020-13). When you use DeerLab in your work, please cite 

> Fábregas Ibáñez, L., Jeschke, G., and Stoll, S.: DeerLab: A comprehensive toolbox for analyzing dipolar EPR spectroscopy data, Magn. Reson. Discuss., https://doi.org/10.5194/mr-2020-13, in review, 2020

Please check back frequently for updated publication information.

### License

The DeerLab toolbox is licensed under the MIT License. The complete toolbox consists of the functions ([functions/](https://github.com/JeschkeLab/DeerLab/tree/master/functions)), documentation source ([docsrc/](https://github.com/JeschkeLab/DeerLab/tree/master/docsrc)), tutorial scripts ([tutorials/](https://github.com/JeschkeLab/DeerLab/tree/master/tutorials)), test suite ([tests/](https://github.com/JeschkeLab/DeerLab/tree/master/tests)), and pipeline scripts ([.github/workflows](https://github.com/JeschkeLab/DeerLab/tree/master/.github/workflows)). See below for exceptions.

DeerLab includes code from the following projects, which have their own licenses:
- [datahash.m](https://www.mathworks.com/matlabcentral/fileexchange/31272-datahash) (Hash-key generator by Jan Simon) [BSD] 
- [fresnelS.m, fresnelC.m](https://www.mathworks.com/matlabcentral/fileexchange/28765-fresnels-and-fresnelc) (Efficient and accurate Fresnel integrals by John D'Errico) [BSD]
- [nlsqbnd.m](https://ch.mathworks.com/matlabcentral/fileexchange/23621-nlsqbnd) (Non-linear least squares solver with box constraints by Alain Barraud) [BSD]
- [golden.m](https://www.mathworks.com/matlabcentral/fileexchange/25919-golden-section-method-algorithm) (Golden Section method algorithm by Katarzyna Zarnowiec) [BSD]
- [jacobianest.m](https://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation) (Adaptive Robust Numerical Differentiation by John D'Errico) [BSD]
- [kde.m](https://ch.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator) (Kernel Density Estimator by Zdravko Botev) [BSD]
- [LevenbergMarquardt.m, jacobiansimple.m](https://ch.mathworks.com/matlabcentral/fileexchange/53449-levenberg-marquardt-toolbox) (Levenberg-Marquardt & Jacobian toolbox by Alexander Dentler)[BSD]
- [fdcoeffF.m](https://faculty.washington.edu/rjl/fdmbook/matlab/fdcoeffF.m) (Fornberg's method for finite difference coefficients)
- [minq.m](http://www.mat.univie.ac.at/~neum/software/minq/) (MINQ8 - General Definite and Bound Constrained Indefinite Quadratic Programming by Waltraud Huyer and Arnold Neumair)

Copyright (c) 2019-2020: Luis Fábregas Ibáñez, Stefan Stoll, Gunnar Jeschke, and [other contributors](https://github.com/JeschkeLab/DeerLab/contributors).
