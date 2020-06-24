.. highlight:: matlab
.. _snlls:

*********************
:mod:`snlls`
*********************

Separable Non-linear Least Squares Solver

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    [pnlin,plin,puq,stats] = snlls(___)
    __ = snlls(y,Amodel,par0,lb,ub,lbl,ubl)
    __ = snlls(y,Amodel,par0,lb,ub,lbl)
    __ = snlls(y,Amodel,par0,lb,ub)
    __ = snlls(y,Amodel,par0,lb)
    __ = snlls(y,Amodel,par0)
    __ = snlls(___,'Name',Values,___)


Parameters
    *   ``y`` - Input data (*N*-element array)
    *   ``Amodel`` - Non-linear function accepting ``W`` non-linear parameters, returns a *NxM*-matrix ``A`` (function handle)
    *   ``par0`` -  Start values of the non-linear parameters (*W*-element array)
    *   ``lb`` -  Lower bounds of the non-linear parameters (*W*-element array)
    *   ``ub`` -  Upper bounds of the non-linear parameters (*W*-element array)
    *   ``lbl`` -  Lower bounds of the linear parameters (*M*-element array)
    *   ``ubl`` -  Upper bounds of the linear parameters (*M*-element array)

Returns
    *  ``pnlin`` - Fitted non-linear parameters
    *  ``plin`` - Fitted linear parameters
    *  ``puq`` - Uncertainty quantification (struct)
    *  ``stats`` - Goodness-of-fit statistics (struct)

-----------------------------


Description
=========================================

This separable non-linear least squares (SNLLS) solver aims to solve the following general problem 


.. code-block:: none 

     [pnlin,plin] = argmin ||Amodel(pnlin)*plin - y||^2
                    subject to  pnlin in [lb,ub]
                                 plin in [lbl,ubl]
 
where the parameter space is composed of a set of non-linear parameters ``pnlin`` and linear parameters ``plin``. If the non-linear function ``Amodel`` yields an ill-conditioned problem, the solver will include a regularization penalty and solve the following problem


.. code-block:: none 

     [pnlin,plin] = argmin ||Amodel(pnlin)*plin - y||^2 + alpha^2*||L*plin||^2
                    subject to  pnlin in [lb,ub]
                                 plin in [lbl,ubl]
                                 
where ``alpha`` and ``L`` are the regularization parameter and operator, respectively. 

-----------------------------                            

.. code-block:: matlab

    [pnlin,plin] = snlls(y,Amodel,par0)

Fits the input data ``y`` via SNNLS with unconstrained non-linear and linear parameters. The start values of the non-linear parameters ``par0`` must be specified. The function returns the fitted non-linear paramter set ``pnlin`` as well as the fitted linear parameter set ``plin``.


-----------------------------


.. code-block:: matlab

    [pnlin,plin] = snlls(y,Amodel,par0,lb,ub,lbl,ubl)
    [pnlin,plin] = snlls(y,Amodel,par0,lb,ub,lbl)
    [pnlin,plin] = snlls(y,Amodel,par0,lb,ub)
    [pnlin,plin] = snlls(y,Amodel,par0,lb)

The boundaries for the non-linear paramters (``lb`` and ``ub``) as well as for the linear parameter (``lbl`` and ``ubl``) can be specified as additional input arguments. If not specified or passed empty, the boundaries are set to infinity (unbounded).

-----------------------------


.. code-block:: matlab

    [pnlin,plin,puq] = snlls(___)

The third output argument contains the uncertainty quantification structure for the full parameter set (non-linear + linear) based on the covariance matrix of the SNNLS problem. In addition to the functionality described in :ref:`cireference`, when requesting the confidence intervals via the ``puq.ci`` field, an additional argument can be passed to request the confidence intervals of the individual linear or non-linear parameter sets.

    *  ``puq.ci(n)`` - Confidence interval of the combined parameter set
    *  ``puq.ci(n,'lin')`` - Confidence interval of the linear parameter set
    *  ``puq.ci(n,'nonlin')`` - Confidence interval of the non-linear parameter set

-----------------------------

.. code-block:: matlab

    [pnlin,plin,puq,stats] = snlls(___)
    
The ``stats`` structure provides several statistical metric which allow judgment on the quality of the fitted ``yfit`` on the experimental data ``y`` and allows comparison between fits. The structure contains the following fields: 

         *   ``.chi2red`` - Reduced `\chi^2` test
         *   ``.R2`` - `R^2` test
         *   ``.RMSD`` - Root-mean squared deviation (RMSD)
         *   ``.AIC`` - Akaike information criterion
         *   ``.AICc`` - Corrected Akaike information criterion
         *   ``.BIC`` - Bayesian information criterion

Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.


.. code-block:: matlab

    ___ = snlls(___,'Name1',Value1,'Name2',Value2,___)


- ``'includePenalty'`` - Include regularization penalty term
    Manually specifies whether a regularization penalty is included (``true``) or not (``false``). If not specified or passed empty, the penalty is automatically included when the matrix returned by ``Amodel`` is ill-conditioned (default). 

    *Default:* []

    *Example:*

		.. code-block:: matlab

			___ = snlls(___,'includePenalty',false)


- ``'alphaOptThreshold'`` - Relative parameter change threshold
    Specifies the relative parameter change threshold for reoptimizing the regularization parameter during the fitting.

    *Default:* 1e-3

    *Example:*

		.. code-block:: matlab

			___ = snlls(___,'alphaOptThreshold',1e-4)


- ``'RegParam'`` - Regularization parameter
    Specifies the selection method employed for the optimization of the regularization parameter (see :ref:`selregparam` for a list). If a value is passed, the regularization parameter will be fixed througout the optimization.

    *Default:* ``'aic'``

    *Example:*

		.. code-block:: matlab

			___ = snlls(___,'RegParam','aic')


- ``'RegType'`` - Regularization functional type
    Type of regularization penalty

        *   ``'tikhonov'`` - Tikhonov regularization
        *   ``'tv'`` - Total variation regularization
        *   ``'huber'`` - Huber regularization

    *Default:* ``tikhonov``

    *Example:*

		.. code-block:: matlab

			___ = snlls(___,'RegType','tv')

- ``'RegOrder'`` - Regularization matrix order
    Order of the regularization operator matrix.

    *Default:* ``2``

    *Example:*

		.. code-block:: matlab

			___ = snlls(___,'RegOrder',0)


- ``'NonLinSolver'`` - Optimization solver for non-linear part 
    Numerical solver employed for fitting non-linear parameters to the data.

        *   ``'lsqnonlin'`` - Non-linear least squares (requires Optimization toolbox)
        *   ``'lmlsqnonlin'`` - Levenberg-Marquardt non-linear least squares (free)

    *Default:* ``'lsqnonlin'`` (Optimization Toolbox installed) or ``'lmlsqnonlin'`` (Optimization Toolbox not installed)

    *Example:*

		.. code-block:: matlab

			___ = snlls(___,'nonLinSolver','lmlsqnonlin')


- ``'LinSolver'`` - Optimization solver for linear part 
    Numerical solver employed for fitting linear parameters to the data.

        *   ``'lsqlin'`` - Linear least squares (requires Optimization toolbox)
        *   ``'minq'`` - Quadratic programming solver [MINQ5](https://www.mat.univie.ac.at/~neum/software/minq/) (free)

    *Default:* ``'lsqlin'`` (Optimization Toolbox installed) or ``'minq'`` (Optimization Toolbox not installed)

    *Example:*

		.. code-block:: matlab

			___ = snlls(___,'LinSolver','lsqlin')
            

- ``'NonLinTolFun'`` -  Tolerance value for the non-linear solver
    Non-linear optimizer function tolerance. The solver stops once the fitting functional evaluation reaches a value lower than this tolerance. Lower values increase the precision of the result, albeit at the cost of longer computation times.

    *Default:* ``1e-5``

    *Example:*

		.. code-block:: matlab

			param = snlls(___,'NonLinTolFun',1e-20)

- ``'NonLinMaxIter'`` - Maximal iterations for the non-linear solver 
    Maximum number of iterations of the non-linear solver. After the solver exceeds this number the optimization will stop.solvers.

    *Default:* ``1e4``

    *Example:*

		.. code-block:: matlab

			param = snlls(___,'nonLinMaxIter',1e10)


- ``'LinTolFun'`` -  Tolerance value for the linear solver
    Linear optimizer function tolerance. The solver stops once the fitting functional evaluation reaches a value lower than this tolerance. Lower values increase the precision of the result, albeit at the cost of longer computation times.

    *Default:* ``1e-5``

    *Example:*

		.. code-block:: matlab

			param = snlls(___,'LinTolFun',1e-20)

- ``'LinMaxIter'`` - Maximal iterations for the linear solver 
    Maximum number of iterations of the linear solver. After the solver exceeds this number the optimization will stop.solvers.

    *Default:* ``1e4``

    *Example:*

		.. code-block:: matlab

			param = snlls(___,'LinMaxIter',1e10)
            
- ``'MultiStart'`` -  Multi-start global optimization
    Number of initial points to be generated for a global search. For each start point, a local minimum is searched, and the solution with the lowest objective function value is selected as the global optimum.

    *Default:* ``1`` (No global optimization)

    *Example:*

		.. code-block:: matlab

			param = snlls(___,'MultiStart',50)

