.. highlight:: matlab
.. _fitmultimodel:


***********************
:mod:`fitmultimodel`
***********************

Dipolar signal fitting with automatic multi-component parametric distance distributions optimization

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    [Pfit,parfit,Puq,paruq,Nopt,fcnal,Peval,stats] = fitmultimodel(___)
    [___] = fitmultimodel(V,K,r,@dd_model,Nmax,metric,lb,ub)
    [___] = fitmultimodel(V,K,r,@dd_model,Nmax,metric,lb)
    [___] = fitmultimodel(V,K,r,@dd_model,Nmax,metric)
    [___] = fitmultimodel(V,@Kmodel,r,@dd_model,Nmax,metric,lb,ub)
    [___] = fitmultimodel(V,@Kmodel,r,@dd_model,Nmax,metric,lb)
    [___] = fitmultimodel(V,@Kmodel,r,@dd_model,Nmax,metric)
    [___] = fitmultimodel({V1,V2,___},{K1,K2,___},r,@dd_model,Nmax,metric,lb,ub)
    [___] = fitmultimodel({V1,V2,___},Kmodel,r,@dd_model,Nmax,metric,lb,ub)
    [___] = fitmultimodel(___,'Name',Value)


Parameters
    *   ``V`` - Input signal (*N*-element array)
    *   ``K`` -  Dipolar kernel (*NxM*-element array)
    *   ``@Kmodel`` -  Kernel model function (function handle)
    *   ``r`` -  Distance Axis (*M*-element array)
    *   ``@dd_model`` -  Parametric model basis function (function handle)
    *   ``Nmax`` - Maximum number of components (scalar)
    *    ``metric`` - Metric for model selection (string)
    *    ``lb`` - Lower bounds for parameters ``{Kpar_lb,ddpar_lb}`` (2-element cell array)
    *    ``ub`` - Upper bounds for parameters ``{Kpar_ub,ddpar_ub}`` (2-element cell array)


Returns
    *  ``Pfit`` - Fitted distance distribution (*M*-element array)
    *  ``parfit`` - Fitted model parameters ``{Kpar,ddpar,Amp}``(3-element cell array)
    *  ``Puq`` - Fitted distribution uncertainty quantification (struct)
    *  ``paramuq`` - Fitted parameters uncertainty quantification (struct)
    *  ``Nopt`` - Optimal number of components (scalar)
    *  ``metrics`` - Evaluated model selection functional (``Nmax``-element array)
    *  ``Peval`` - Fitted distance distributions for each multi-component model (``Nmax``-element cell array)

-----------------------------


Description
=========================================

This function fits a multi-component distance distribution (of maximal ``Nmax`` components) via the problem  

.. code-block:: none 

     [Amp,Kpar,ddpar] = argmin || sum_i Amp_i*Kmodel(Kpar)*dd_model(ddpar_i) - V||^2
                        subject to   Kpar  in [Kpar_lb,Kpar_ub]
                                     ddpar in [ddpar_lb,ddpar_ub]
                                     A_i   in [0,inf] 

where ``Kpar`` (kernel parameters) and ``ddpar`` (distribution parameters) are fitted as the non-linear parameters and ``Amp`` (components amplitudes) as the linear parameters of a separable non-linear least-squares (SNLLS) problem. The optimization is repeated for all multi-component models up to ``Nmax``, and the optimal model is chosen according to a selection metric specified by ``metric``.

-----------------------------

.. code-block:: matlab

        Pfit = fitmultimodel(V,K,r,@dd_model,Nmax,metric)
        Pfit = fitmultimodel(V,K,r,@dd_model,Nmax,metric,{[],ddpar_lb},{[],ddpar_ub})


Fits a distance distribution ``Pfit`` to the dipolar signal ``V`` using a multi-component parametric model and a fixed dipolar kernel ``K``. The multi-component model is constructed from the basis function provided as ``@dd_model``. The function chooses the optimal number of components distributions up to a maximum number given by ``Nmax`` by means of the metric specified as the input ``metric``. The accepted inputs are:

	*   ``'aic'`` - Akaike information criterion
	*   ``'aicc'`` - Corrected Akaike information criterion
	*   ``'bic'`` - Bayesian information criterion
	*   ``'rmsd'`` - Root mean squared deviation

The lower and upper bounds of the parameters accepted by ``@dd_model`` can be specified as the ``ddpar_lb`` and ``ddpar_ub`` arguments. If not specified, they are automatically taken from the model defition. 

-----------------------------

.. code-block:: matlab

        Pfit = fitmultimodel(V,@Kmodel,r,@dd_model,Nmax,metric)
        Pfit = fitmultimodel(V,@Kmodel,r,@dd_model,Nmax,metric,{Kpar_lb,[]},{Kpar_ub,[]})
        Pfit = fitmultimodel(V,@Kmodel,r,@dd_model,Nmax,metric,{Kpar_lb,ddpar_lb},{Kpar_ub,ddpar_ub})

If the kernel depends on other parameters (e.g. modulation depth, background parameters,...) these can be fitted along the multi-component distribution by passing a custom dipolar kernel model as a function handle ``Kmodel``. This must be a function that accepts an array ``Kpar`` of parameters and returns a *NxM*-element dipolar kernel matrix. For example: 

 .. code-block:: matlab

        Kmodel = @(Kpar) K4pdeer(Kpar,t,r);
        Pfit = fitmultimodel(V,Kmodel,r,@dd_model,Nmax,'aicc')
        
        function K = K4pdeer(Kpar,t,r)
            lambda = Kpar(1);
            conc = Kpar(2);
            K = dipolarkernel(t,r,lambda,bg_hom3d(t,conc,lambda));
        end        

The lower and upper bounds of the parameters accepted by ``@Kmodel`` can be specified as the ``Kpar_lb`` and ``Kpar_ub`` arguments. If not specified, they are considered to be unconstrained. The fitted kernel parameters are returned along the distribution parameters in the ``parfit`` output.

-----------------------------


.. code-block:: matlab

    P = fitmultimodel({V1,V2,___},{K1,K2,___},r,@dd_model,Nmax,metric)
    P = fitmultimodel({V1,V2,___},Kmodel,r,@dd_model,Nmax,metric)


Passing multiple signals results in global fitting of a single multi-component distance distribution. The multiple signals are passed as a cell array of arrays of sizes *N1*, *N2*,... and if a fixed kernel is used, then a cell array of kernel matrices with sizes *N1xM*, *N2xM*, ... must be passed as well. If instead a kernel model is used, a single function handle ``Kmodel`` must be specified, whose function returns a cell array of kernels.  For example:


 .. code-block:: matlab

        Kmodel = @(Kpar) K4pdeer(Kpar,t1,t2,r);
        Pfit = fitmultimodel({V1,V2,___},Kmodel,r,@dd_model,Nmax,'aicc')
        
        function K = K4pdeer(Kpar,t1,t2,r)
            lambda1 = Kpar(1);
            conc1 = Kpar(2);
            lambda2 = Kpar(1);
            conc2 = Kpar(2);
            K{1} = dipolarkernel(t1,r,lambda1,bg_hom3d(t1,conc1,lambda1));
            K{2} = dipolarkernel(t2,r,lambda2,bg_hom3d(t2,conc2,lambda2));
        end    

-----------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.


.. code-block:: matlab

    ___ = fitmultimodel(___,'Property1',Value1,'Property2',Value2,___)

- ``'GlobalWeights'`` - Global analysis weights
    Array of weighting coefficients for the individual datasets in global fitting. If not specified, the global fit weights are automatically computed according to their contribution to ill-posedness. The same number of weights as number of input signals is required. Weight values do not need to be normalized.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			___ = fitmultimodel({V1,V2,V3},___,'GlobalWeights',[0.1 0.6 0.3]])
            

- ``'normP'`` -  Renormalization of the distance distribution
    This enables/disables the re-normalization of the fitted distance distribution such that ``trapz(r,P)==1``. 

    *Default:* ``true``

    *Example:*

		.. code-block:: matlab

			P = fitregmodel(___,'normP',false)


- See :ref:`snlls` for a detailed list of other property-value pairs accepted by the function.