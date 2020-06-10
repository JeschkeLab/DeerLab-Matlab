.. highlight:: matlab
.. _dd_shellvoidsphere:


***************************
:mod:`dd_shellvoidsphere`
***************************

Particles distributed on a sphere inside a spherical shell separated by a void 

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_shellvoidsphere()
        P = dd_shellvoidsphere(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)

-----------------------------

Model
=========================================

.. image:: ../images/model_scheme_dd_shellvoidsphere.png
   :width: 25%

:math:`P(r) = ((R_3^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_3) - (R_2^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_2) )/(R_3^3 - R_2^3)`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + d`

:math:`R_3 = R + d + w`

================ ============== ========= ======== ========= ===================================
 Variable         Symbol         Default   Lower    Upper       Description
================ ============== ========= ======== ========= ===================================
``param(1)``     :math:`R`       1.5       0.1        20        Sphere radius
``param(2)``     :math:`w`       1.0       0.1        20        Shell thickness
``param(2)``     :math:`d`       0.5       0.1        20        Shell-Sphere separation
================ ============== ========= ======== ========= ===================================


Example using default parameters:

.. image:: ../images/model_dd_shellvoidsphere.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_shellvoidsphere()

Returns an ``info`` structure containing the information of the model parameters and boundaries.

* ``info(n).Index`` -  Index of the parameter in the ``param`` array.
* ``info(n).Parameter`` -  Description of the n-th parameter.
* ``info(n).Lower`` -  Lower bound of the n-th parameter.
* ``info(n).Upper`` -  Upper bound of the n-th parameter.
* ``info(n).Start`` -  Start value of the n-th parameter.

-----------------------------


.. code-block:: matlab

    P = dd_shellvoidsphere(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

References
=========================================

[1] D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63.
`DOI:  10.1016/j.jmr.2013.01.007 <http://doi.org/10.1016/j.jmr.2013.01.007>`_