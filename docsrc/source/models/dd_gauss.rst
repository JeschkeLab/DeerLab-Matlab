.. highlight:: matlab
.. _dd_gauss:


***********************
:mod:`dd_gauss`
***********************

Gaussian distribution

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_gauss()
        P = dd_gauss(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (table)

-----------------------------

Model
=========================================

:math:`P(r) = \sqrt{\frac{2}{\pi}}\frac{1}{\sigma}\exp\left(-\frac{(r-\left<r\right>)^2}{\sigma^2}\right)`

with :math:`\sigma = \mathrm{FWHM}/\sqrt{2ln(2)}`

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`r0`                 3.5         1.0              20         center distance (nm)
``param(2)``   :math:`\mathrm{FWHM}`      0.5         0.2              5          FWHM (nm)
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_gauss.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_gauss()

Returns an ``info`` table containing the information of the model parameters and boundaries.

The table contents can be accessed as follows:
* ``info.Index`` -  Indices of the parameters in the ``param`` array
* ``info.Parameter`` -  Names of the model parameters
* ``info.Lower`` - Lower bounds for the parameters
* ``info.Upper`` - Upper bounds for the parameters
* ``info.Start`` - Start values for optimization

-----------------------------


.. code-block:: matlab

    P = dd_gauss(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

