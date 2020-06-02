.. highlight:: matlab
.. _bootan:

***********************
:mod:`bootan`
***********************

Bootstrap analysis for uncertainty estimation

------------------------


Syntax
=========================================

.. code-block:: matlab

    [bootci,stats] = bootan(fcn,V,Vfit)
    [bootci,stats] = bootan(fcn,V,Vfit,Nsamples)
    [bootci,stats] = bootan(___,'Property',Value)

Parameters
    *   ``fcn`` - Function to analyze (function handle)
    *   ``V`` - Experimental signal (*N*-element array)
    *   ``Vfit`` - Fitted signal (*N*-element array)
    *   ``Nsamples`` - Number of bootstrap samples (scalar)

Returns
    *   ``bootci`` - Bootstrapped confidence intervals (struct or cell array)

------------------------


Description
=========================================

.. code-block:: matlab

    [bootci,stats] = bootan(fcn,V,Vfit)

Performs a uncertainty analysis of the output variables of the function ``fcn`` from 1000 bootstrap samples. The output argument of ``fcn`` is evaluated for all level-combinations of the factors. The function to be analyzed must be a function handle accepting the ``V`` experimental signal as input. Example:

.. code-block:: matlab

    [bootci,stats] = bootan(@(V)myfcn(p,varargin),V,Vfit)

    function [Pfit1,Pfit2] = myfcn(V,varargin)
       Pfit1 = fitparamodel(V,@dd_gauss,r,K)
       Pfit2 = fitparamodel(V,@dd_randcoil,r,K)
    end


The output ``bootci`` is a single (in case of a single evaluated variable) or a cell array (in case of multiple evaluated variables) of confidence intervals structures (see :ref:`cireference`).

------------------------

.. code-block:: matlab

    stats = bootan(fcn,V,Vfit,Nsamples)


The number of bootstrap samples can be specified in ``Nsamples``. The quality of bootstrapping results improve with the number of boostrap samples evaluated. 



------------------------


.. code-block:: matlab

    stats = bootan(fcn,{V1,V2,___},{Vfit1,Vfit2,___},Nsamples)


If the evaluated function ``fcn`` requries multiple signals ``{V1,V2,___}`` as input, these can be specified aloong the same number of fitted signals ``{Vfit1,Vfit2,___}``. 


------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.



.. code-block:: matlab

    stats = bootan(___,'Property1',Value1,'Property2',Value2,___)

- ``'Verbose'`` - Display progress information
    Specifies whether to print the progress of the bootstrap analysis on the command window.

    *Default:* ``false``

    *Example:*

		.. code-block:: matlab

			stats = bootan(___,'Verbose',true)


- ``'Resampling'`` - Re-sampling method
    Specifies the method employed for re-sampling new bootstrap samples.

        *   ``'gaussian'`` - Sample noise from a Gaussian distribution
        *   ``'residual'`` - Sample noise from the fit residuals

    *Default:* ``gaussian``

    *Example:*

		.. code-block:: matlab

			stats = bootan(___,'Resampling',residual)

