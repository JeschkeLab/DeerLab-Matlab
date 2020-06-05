.. highlight:: matlab
.. _dd_uniform:


***********************
:mod:`dd_uniform`
***********************

Uniform distribution

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_uniform()
        P = dd_uniform(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (table)

-----------------------------

Model
=========================================


This provides a simple uniform distribution.

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`r_\mathrm{L}`         2.5       0.1              6           left edge
``param(2)``   :math:`r_\mathrm{R}`         3.0       0.2              20          right edge
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_uniform.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_uniform()

Returns an ``info`` table containing the information of the model parameters and boundaries.

The table contents can be accessed as follows:
* ``info.Index`` -  Indices of the parameters in the ``param`` array
* ``info.Parameter`` -  Names of the model parameters
* ``info.Lower`` - Lower bounds for the parameters
* ``info.Upper`` - Upper bounds for the parameters
* ``info.Start`` - Start values for optimization


-----------------------------


.. code-block:: matlab

    P = dd_uniform(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

