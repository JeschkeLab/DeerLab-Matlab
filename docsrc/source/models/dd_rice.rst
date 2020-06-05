.. highlight:: matlab
.. _dd_rice:


***********************
:mod:`dd_rice`
***********************

3D-Rice distribution

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_rice()
        P = dd_rice(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (table)

-----------------------------

Model
=========================================

:math:`P(r) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=3` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`. This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`\nu`                3.5     1.0              10         Center
``param(2)``   :math:`\sigma`             0.7     0.1              5          Width
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_rice.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_rice()

Returns an ``info`` table containing the information of the model parameters and boundaries.

The table contents can be accessed as follows:
* ``info.Index`` -  Indices of the parameters in the ``param`` array
* ``info.Parameter`` -  Names of the model parameters
* ``info.Lower`` - Lower bounds for the parameters
* ``info.Upper`` - Upper bounds for the parameters
* ``info.Start`` - Start values for optimization

-----------------------------


.. code-block:: matlab

    P = dd_rice(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

