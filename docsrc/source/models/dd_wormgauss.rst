.. highlight:: matlab
.. _dd_wormgauss:

***********************
:mod:`dd_wormgauss`
***********************

Worm-like chain model near the rigid limit with Gaussian convolution

Syntax
=========================================

.. code-block:: matlab

        info = dd_wormgauss()
        P = dd_wormgauss(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (table)

============== =============== ======== ======== ======== ===============================
 Variable       Symbol         Default   Lower   Upper       Description
============== =============== ======== ======== ======== ===============================
``param(1)``   :math:`L`       3.7      1.5       10        Contour length
``param(2)``   :math:`L_p`     10       2         100       Persistence length
``param(3)``   :math:`\sigma`  0.2      0.01      2         Gaussian standard deviation
============== =============== ======== ======== ======== ===============================

Example using default parameters:

.. image:: ../images/model_dd_wormgauss.png
   :width: 650px


Description
=========================================

.. code-block:: matlab

        info = dd_wormgauss()

Returns an ``info`` table containing the information of the model parameters and boundaries.

The table contents can be accessed as follows:
* ``info.Index`` -  Indices of the parameters in the ``param`` array
* ``info.Parameter`` -  Names of the model parameters
* ``info.Lower`` - Lower bounds for the parameters
* ``info.Upper`` - Upper bounds for the parameters
* ``info.Start`` - Start values for optimization

.. code-block:: matlab

    P = dd_wormgauss(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.



References
=========================================

[1] J. Wilhelm, E. Frey, Phys. Rev. Lett. 77(12), 2581-2584 (1996)
DOI:  `10.1103/PhysRevLett.77.2581 <https://doi.org/10.1103/PhysRevLett.77.2581>`_