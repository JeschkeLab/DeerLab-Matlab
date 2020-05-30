.. highlight:: matlab
.. _fitparamodel:

*********************
:mod:`fitparamodel`
*********************

Fits a time- or distance-domain parametric model to one (or several) signals.

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    [param,fit,paramci,fitci,stats] = fitparamodel(V,@model,t)
    [param,fit,paramci,fitci,stats] = fitparamodel(V,@model,t,param0)
    [param,fit,paramci,fitci,stats] = fitparamodel({V1,V2,___},@model,{t1,t2,___},param0)
    [param,fit,paramci,fitci,stats] = fitparamodel(V,@model,r,K)
    [param,fit,paramci,fitci,stats] = fitparamodel(V,@model,r,K,param0)
    [param,fit,paramci,fitci,stats] = fitparamodel({V1,V2,___},@model,r,{K1,K2,___},param0)
    [param,fit,paramci,fitci,stats] = fitparamodel(___,'Property',Value)


Parameters
    *   ``V`` - Input signal (*N*-element array)
    *   ``model`` - Parametric model (function handle)
    *   ``t`` -  Model time axis (*N*-element array)
    *   ``r`` -  Model distance axis (*M*-element array)
    *   ``K`` -  Dipolar kernel (*NxM*-element array)
    *   ``param0`` -  Model parameter initial guess (*W*-array)
Returns
    *  ``param`` - Fitted model parameters (*W*-array)
    *  ``fit`` - Parametric model fit (*N*-element array)
    *  ``paramci`` - Fit confidence intervals (*Wx2*-element array)
    *  ``fitci`` - Model fit confidence bands (*Nx2*-element array)
    *  ``stats`` - Goodness-of-fit statistics (structure)


-----------------------------


Description
=========================================

.. code-block:: matlab

    [param,fit,paramci,fitci] = fitparamodel(V,@model,t)
    [param,fit,paramci,fitci] = fitparamodel(V,@model,t,param0)

Fits the **time-domain** parametric model ``@model`` to the input signal ``V`` on a time axis ``t``. User-defined initial guess values can be passed as an additional argument, if not they are automatically determined from the model. If the model is a user-defined function handle, the function will require ``param0`` to be passed.

The fitted parameters as well as the corresponding confidence intervals structures (see :ref:`cireference`) are returned as the ``param`` and ``paramci`` outputs, respectively. The fitted model and its confidence bands are returned as the ``fit`` and ``fitci`` outputs, respectively.


-----------------------------


.. code-block:: matlab

    [param,fit,paramci,fitci] = fitparamodel(V,@model,r,K)
    [param,fit,paramci,fitci] = fitparamodel(V,@model,r,K,param0)

Fits the **distance-domain** parametric model ``@model`` to the input signal ``V`` on a distance axis ``r``. The dipolar kernel ``K`` is required as in input for distance-domain fitting. User-defined initial guess values can be passed as an additional argument, if not they are automatically determined from the model. If the model is a user-defined function handle, the function will require ``param0`` to be passed.

-----------------------------


.. code-block:: matlab

    param = fitparamodel({V1,V2,___},@model,r,{K1,K2,___})
    param = fitparamodel({V1,V2,___},@model,r,{K1,K2,___},param0)

Passing multiple signals/kernels enables **distance-domain global fitting** of the parametric model to a single distribution. The global fit weights are automatically computed according to their contribution to ill-posedness. The multiple signals are passed as a cell array of arrays of sizes *N1*, *N2*,... and a cell array of Kernel matrices with sizes *N1xM*, *N2xM*, ... must be passed as well.

-----------------------------


.. code-block:: matlab

    param = fitparamodel({V1,V2,V3},@model,{t1,t2,t3})
    param = fitparamodel({V1,V2,V3},@model,{t1,t2,t3},param0)

Similarly, **time-domain global fitting** can be used when passing a time-domain ``@model`` and the model time axes ``{t1,t2,___}`` of the corresponding signals.

-----------------------------


User-defined parametric models must have the following function definition structure:

.. code-block:: matlab

    Vfit = model(t,param)
    Pfit = model(r,param)
	
where the ``r`` and ``t`` depend on whether the parametric model is a distance or time-domain model, respectively. Additionally the parametric model can accept a third input argument ``idx`` as follows

.. code-block:: matlab

    Vfit = model(t,param,idx)
    Pfit = model(r,param,idx)

By doing so, ``fitparamodel`` will automatically pass the index ``idx = (1,2,...,N)`` of the input signal cell array  
``{S1,S2,...,SN}`` being currently processed. This allows for implementation different routines in the parametric model for different signals during global fitting. 


-----------------------------

.. code-block:: matlab

    [param,fit,paramci,fitci,stats] = fitparamodel(___)

The ``stats`` structure provides several statistical metric which allow judgment on the quality of the fitted ``Vfit`` on the experimental data ``V`` and allows comparison between fits. The structure contains the following fields: 

         *   ``.chi2red`` - Reduced `\chi^2` test
         *   ``.R2`` - `R^2` test
         *   ``.RMSD`` - Root-mean squared deviation (RMSD)
         *   ``.AIC`` - Akaike information criterion
         *   ``.AICc`` - Corrected Akaike information criterion
         *   ``.BIC`` - Bayesian information criterion

-----------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.


.. code-block:: matlab

    param = fitparamodel(___,'Property1',Value1,'Property2',Value2,___)


- ``'Upper'`` - Parameters upper bound constraints
    Array of upper bounds for the model parameters.

    *Default:* unbounded or automatically set

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(___,'Upper',[1 100])

- ``'Lower'`` - Parameters lower bound constraints
    Array of lower bounds for the model parameters.

    *Default:* unbounded or automatically set

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(___,'Lower',[0 3])

- ``'Solver'`` - Optimization solver
    Numerical solver employed for fitting the model to the data.

        *   ``'lsqnonlin'`` - Non-linear least squares (requires Optimization toolbox)
        *   ``'lmlsqnonlin'`` - Levenberg-Marquardt non-linear least squares (free)
        *   ``'nlsqbnd'`` - Non-linear least squares (free, Windows OS only)

    *Default:* ``'lsqnonlin'`` (Optimization Toolbox installed) or ``'lmlsqnonlin'`` (Optimization Toolbox not installed)

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(___,'Solver','lmlsqnonlin')

- ``'Algorithm'`` - Numerical solver algorithm
    Algorithm to be used by the solvers (see ``lsqnonlin`` MATLAB documentation)

    *Default:* see MATLAB documentation

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(___,'Algorithm','trust-region-reflective')

- ``'GlobalWeights'`` - Global analysis weights
    Array of weighting coefficients for the individual signals in global fitting. If not specified, the global fit weights are automatically computed according to their contribution to ill-posedness. The same number of weights as number of input signals is required. Weight values do not need to be normalized.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			param = fitparamodel({S1,S2,S3},@dd_gauss,r,{K1,K2,K3},'GlobalWeights',[0.1 0.6 0.3]])

- ``'TolFun'`` -  Optimizer tolerance value
    Optimizer function tolerance. The solver stops once the fitting functional evaluation reaches a value lower than this tolerance. Lower values increase the precision of the result, albeit at the cost of longer computation times.

    *Default:* ``1e-9``

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(___,'TolFun',1e-20)

- ``'MaxIter'`` - Maximal solver iterations
    Maximum number of iterations of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'``  and ``'lsqnonneg'`` solvers.

    *Default:* ``2e7``

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(___,'MaxIter',1e10)

- ``'MaxFunEval'`` -  Maximal solver function evaluations
    Maximum number of function evaluation of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'``  and ``'lsqnonneg'`` solvers.

    *Default:* ``2e7``

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(___,'MaxFunEval',1e10)

- ``'Rescale'`` -  Rescaling of fitted dipolar signal
    This enables/disables the automatic optimization of the dipolar signal scale. If enabled (``true``) the experimental dipolar signal does not need to fulfill ``V(t=0)=1``, if disabled (``false``) it needs to be fulfilled.

    *Default:* ``true``

    *Example:*

		.. code-block:: matlab

			V = correctscale(V,t);
			param = fitparamodel(___,'Rescale',false)

- ``'MultiStart'`` -  Multi-start global optimization
    Number of initial points to be generated for a global search. For each start point, a local minimum is searched, and the solution with the lowest objective function value is selected as the global optimum.

    *Default:* ``1`` (No global optimization)

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(___,'MultiStart',50)

- ``'Verbose'`` -  Information display
    Set the level of detail display for the solvers:

        *   ``'off'`` - No information displayed
        *   ``'final'`` - Display solver exit message
        *   ``'iter-detailed'`` - display state of solver at each iteration


    *Default:* ``'off'``

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(___,'Verbose','iter-detailed')