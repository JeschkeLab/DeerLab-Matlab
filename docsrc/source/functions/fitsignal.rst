.. highlight:: matlab
.. _fitsignal:

***********************
:mod:`fitsignal`
***********************

Fit a full model to a dipolar time-domain trace

------------------------


Syntax
=========================================

.. code-block:: matlab

    fitsignal(V,t,r,dd,bg,ex,par0)
    [Vfit,Pfit,Bfit,parfit,modfitci,parci,stats] = fitsignal(V,t,r,dd,bg,ex,par0,lb,ub)
    __ = fitsignal(V,t,r,dd,bg,ex,par0)
    __ = fitsignal(V,t,r,dd,bg,ex)
    __ = fitsignal(V,t,r,dd,bg)
    __ = fitsignal(V,t,r,dd)
    __ = fitsignal(V,t,r)
    __ = fitsignal({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},{ex1,ex2,__},par0,lb,ub)
    __ = fitsignal({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},{ex1,ex2,__},par0)
    __ = fitsignal({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},{ex1,ex2,__})
    __ = fitsignal({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},ex)
    __ = fitsignal({V1,V2,__},{t1,t2,__},r,dd)
    __ = fitsignal({V1,V2,__},{t1,t2,__},r)
    __ = fitsignal(___,'Property',Values,___)

Parameters
    *   ``V`` -- Time-domain signal to fit (*N*-element array)
    *   ``t`` -- Time axis, in microseconds (*N*-element array)
    *   ``r`` -- Distance axis, in nanometers (*M*-element array)
    *   ``dd`` -- Distance distribution model, can be...

                 * ``@dd_model`` function handle for parametric distribution model
                 * ``'P'`` to indicate a non-parametric distribution
                 * ``'none'`` to indicate no distribution, i.e. only background


    *   ``bg`` -- Background model, can be...

                 -- ``@bg_model`` function handle for parametric background model
                 -- ``'none'`` to indicate no background decay

    *   ``ex`` - Experiment model, can be...

                 -- ``@ex_model`` function handle for parametric experiment model
                 -- ``'none'`` to indicate simple dipolar oscillation (mod.depth = 1)\
                
    *   ``par0`` -- Starting parameters ``{par0_dd,par0_bg,par0_ex}`` (3-element cell array)
    *   ``lb`` -- Lower bounds for parameters ``{lb_dd,lb_bg,lb_ex}`` (3-element cell array)
    *   ``ub`` -- Upper bounds for parameters ``{lb_dd,lb_bg,lb_ex}`` (3-element cell array)


Returns
    *   ``Vfit`` -- Fitted time-domain signal (*N*-element array)
    *   ``Pfit`` -- Fitted distance-domain signal (*M*-element array)
    *   ``Bfit`` -- Fitted background decay (*N*-element array)
    *   ``parfit`` - Structure with fitted parameters (struct)

                 * ``.dd`` -- Fitted parameters for distance distribution model
                 * ``.bg`` -- Fitted parameters for background model
                 * ``.ex`` -- Fitted parameters for experiment model

    *   ``modfituq`` -- Structure with uncertainty quantification of the fitted...
                 * ``.Pfit`` -- ... distance distribution
                 * ``.Bfit`` -- ... background(s)
                 * ``.Vfit`` -- ... signal(s)
                 
    *   ``paruq`` -- Structure with uncertainty quantification of the fitted...
                 * ``.dd`` -- ...distance distribution model parameters
                 * ``.bg`` -- ...background(s) model parameters
                 * ``.ex`` -- ...experiment(s) model parameters
                 
    *   ``stats`` -- Goodness of fit statistical estimators (struct)

------------------------


Description
=========================================

.. code-block:: matlab

    __ = fitsignal(V,t,r,dd,bg,ex)
    __ = fitsignal(V,t,r,dd,bg)
    __ = fitsignal(V,t,r,dd)
    __ = fitsignal(V,t,r)

Fits a full time-domain model of the dipolar signal constructed from the distance distribution model ``dd``, background model ``bg`` and experiment model ``ex`` to the experimental data ``V``, defined on a time-axis ``t``. The distance distribution is fitted on the specified distance axis ``r``. If some models are not specified, the defaults are used: ``'P'`` for the ``dd`` model, ``@bg_hom3d`` for the ``bg`` model, and ``@ex_4pdeer`` for the experiment model.

The fitted dipolar signal ``Vfit``, fitted distribution ``Pfit`` and fitted background ``Bfit`` are returned as the first outputs and their corresponding uncertainty quantification structures as ``modfituq``. The corresponding model parameters are returned in the ``parfit`` structure and their corresponding uncertainty quantification structures (see :ref:`cireference`) in the ``parci`` structure. ``stats`` contains information on the quality of the fit.

.. code-block:: matlab

    fitsignal(V,t,r,dd,bg,ex)

If the function is called without outputs, the function plots the fit results, prints a summary of the fit results, and lists all parameters and their confidence intervals. 

Examples:

.. code-block:: matlab

    fitsignal(V,t,r,@dd_gauss,@bg_hom3d,@ex_4pdeer)  % Fit a 4pDEER signal with homogenous 3D background with Gaussian distribution
    fitsignal(V,t,r,'P',@bg_hom3d,@ex_5pdeer)        % Fit a 5pDEER signal with exponential background and Tikhonov regularization
    fitsignal(V,t,r,'none',@bg_strexp,@ex_4pdeer)    % Fit a 4pDEER stretched exponential background (no foreground)
    fitsignal(V,t,r,@dd_rice,'none','none')          % Fit a dipolar evolution function with Rician distribution
    fitsignal(V,t,r,@dd_gauss2,'none',@ex_4pdeer)    % Fit a 4pDEER form factor (no background) with bimodal Gaussian distribution

------------------------

    __ = fitsignal(V,t,r,dd,bg,ex,par0)

The starting values of the parameter search can be specified via ``par0``. It must be 3-element cell arrays of the form ``{par0_dd,par0_bg,par0_ex}``, where the elements are arrays that give the initial values for the distance distribution parameters, background parameters, and experiment parameters, respectively. If not specified or passed empty, these values are automatically taken from the info structures of the parametric models.


------------------------

    __ = fitsignal(V,t,r,dd,bg,ex,par0,lb,ub)
    __ = fitsignal(V,t,r,dd,bg,ex,par0,[],ub)
    __ = fitsignal(V,t,r,dd,bg,ex,par0,lb)

The lower/upper bounds for the parameter search range can be specified via ``lb`` and ``ub``. These inputs must be 3-element cell arrays of the form ``{ub_dd,ub_bg,ub_ex}`` and ``{lb_dd,lb_bg,lb_ex}``, where the elements are arrays that give the upper/lower bounds for the distance distribution parameters, background parameters, and experiment parameters, respectively. If not specified or passed empty, the boundaries are automatically taken from the info structures of the parametric models.

------------------------

.. code-block:: matlab

    __ = fitsignal({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},{ex1,ex2,__})



Multiple dipolar signals ``{V1,V2,__}`` can be globally fitted to a global distance distribution specified by the model ``dd``. For each signal passed, an experiment and background model can be specified for each signal. The corresponding time-axes ``{t1,t2,__}`` must be provided for all signals respectively.

.. code-block:: matlab

    __ = fitsignal({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},ex)
    __ = fitsignal({V1,V2,__},{t1,t2,__},r,dd,bg,{ex1,ex2,__})
    __ = fitsignal({V1,V2,__},{t1,t2,__},r,dd,bg,ex)
    __ = fitsignal({V1,V2,__},{t1,t2,__},r)

If only one background model ``dd`` or experiment model ``ex`` are specified, that single model is applied for all input signals. If not models are specified, the default models mentioned above are used. 


Examples:

.. code-block:: matlab

    fitsignal({V1,V2},{t1,t2},r,@dd_gauss,@bg_hom3d,{@ex_4pdeer,@ex_4pdeer})  % Fit a Gaussian distribution to a 4pDEER and a 5pDEER signal globally
    fitsignal({V1,V2},{t1,t2},r,'P',@bg_hom3d,@ex_5pdeer)          % Fit a Tikhonov regularized distribution to two different 4pDEER signals

------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the name-value pairs can be passed in any order after the required input arguments have been passed.


.. code-block:: matlab

    fitsignal(___,'Property1',Value1,'Property2',Value2,___)


- ``'RegType'`` - Regularization functional type
    Specifies the type of regularization to be used to fit non-parametric distributions

        *   ``'tikh'`` -- Tikhonov regularization
        *   ``'tv'`` -- Total variation regularization
        *   ``'huber'`` --  Huber regularization

    *Default:* ``tikh``

    *Example:*

		.. code-block:: matlab

			fitsignal(___,'RegType','tv')


- ``'RegParam'`` - Regularization parameter selection
    Specifies the method for the selection of the optimal regularization parameter (``'aic'``, ``'bic'``,...). See ``selregparam`` for more details. The regularization parameter can be manually fixed by passing its value.

    *Default:* ``'aic'``

    *Example:*

		.. code-block:: matlab

			fitsignal(___,'RegParam','bic')
			fitsignal(___,'RegParam',0.2)


- ``'GlobalWeights'`` - Global analysis weights
    Array of weighting coefficients for the individual signals in global fitting. If not specified, the global fit weights are automatically computed according to their contribution to ill-posedness. The same number of weights as number of input signals is required. Weight values do not need to be normalized.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			param = fitparamodel({S1,S2,S3},@dd_gauss,r,{K1,K2,K3},'GlobalWeights',[0.1 0.6 0.3]])

- ``'MultiStart'`` -  Multi-start global optimization
    Number of initial points to be generated for a global search. For each start point, a local minimum is searched, and the solution with the lowest objective function value is selected as the global optimum.

    *Default:* ``1`` (No global optimization)

    *Example:*

		.. code-block:: matlab

			param = fitsignal(___,'MultiStart',50)


- ``'Rescale'`` -  Rescaling of fitted dipolar signal
    This enables/disables the automatic optimization of the dipolar signal scale. If enabled (``true``) the experimental dipolar signal does not need to fulfill ``V(t=0) = 1``, if disabled (``false``) it needs to be fulfilled.

    *Default:* ``true``

    *Example:*

		.. code-block:: matlab

			V = correctscale(V,t);
			fitsignal(___,'Rescale',false)


- ``'normP'`` -  Renormalization of the distance distribution
    This enables/disables the re-normalization of the fitted distance distribution such that ``sum(Pfit)*dr = 1``. 

    *Default:* ``true``

    *Example:*

		.. code-block:: matlab

			fitsignal(___,'normP',false)
            
            
- ``'Display'`` -  Plot and print results
    This enables/disables plotting and printing of results.

    *Default:* ``true`` if no outputs requested, otherwise ``false`` 

    *Example:*

		.. code-block:: matlab

			fitsignal(___,'Display',true)