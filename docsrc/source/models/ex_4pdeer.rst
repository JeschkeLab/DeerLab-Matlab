.. highlight:: matlab
.. _ex_4pdeer:


***********************
:mod:`ex_4pdeer`
***********************

4-pulse DEER experiment 

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = ex_4pdeer()
        pathways = ex_4pdeer(param)

Parameters
    *   ``param`` - Model parameters (array)
Returns
    *   ``pathways`` - Dipolar pathways (array)
    *   ``info`` - Model information (table)


-----------------------------

Model
=========================================

.. image:: ../images/model_scheme_ex_4pdeer.png
   :width: 550px

This experiment model has one modulated pathway and an unmodulated contribution. The kernel is 

.. math::
   K(t,r) =
   [1-\lambda + \lambda K_0(t-T_0^{(1)},r)]B(t-T_0^{(1)},\lambda)

where :math:`T_0^{(1)}=0` is the refocusing time of the modulated dipolar pathway.


============== ================ ============ ============ ============ ================================================
 Variable        Symbol           Default       Lower        Upper                Description
============== ================ ============ ============ ============ ================================================
``param(1)``   :math:`\lambda`     0.3           0            1          modulated pathway amplitude (modulation depth)
============== ================ ============ ============ ============ ================================================


Example of a simulated signal using default parameters:

.. image:: ../images/model_ex_4pdeer.png
   :width: 550px

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = ex_4pdeer()

Returns an ``info`` table containing the information of the model parameters and boundaries.

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    pathways = ex_4pdeer(param)

Generates the dipolar pathways matrix ``pathways`` from the model parameters ``param``. 


