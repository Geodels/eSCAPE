##################
API Documentation
##################

**eSCAPE** is a long-term surface evolution model built to simulate landscape development, sediment transport and
sedimentary basins formation from upstream regions down to marine environments.

Below you will find the *nitty-gritty components* of the Python code...
The following API documentation provides a nearly complete description of the different functions
that compose the code... |:v:|

.. warning::
  It is worth mentioning that a set of **fortran** functions are also part of **eSCAPE** code but are not described
  in this API...



Model class
-----------

.. automodule:: eSCAPE
    :members:

Input/Output
------------

Declaration of input/output functions to

* Parse the YAML input file
* Write HDF5 output of computed surface parameters

.. automodule:: eSCAPE.tools

Input functions
^^^^^^^^^^^^^^^

.. automodule:: eSCAPE.tools.inputparser
    :members:

Output functions
^^^^^^^^^^^^^^^^

.. automodule:: eSCAPE.tools.outmesh
    :members:


Mesher
------------

Set of functions used to build `PETSC` DMPlex mesh and update forcing conditions such as rain, sealevel, tectonic.

.. automodule:: eSCAPE.mesher
    :members:

Unstructured Mesh
^^^^^^^^^^^^^^^^^^

.. automodule:: eSCAPE.mesher.unstructuredmesh
    :members:

Landscape evolution
--------------------

Main functions used to compute surface erosion/deposition. This encompasses both river and hillslope processes.

.. automodule:: eSCAPE.flow
    :members:

Surface processes
^^^^^^^^^^^^^^^^^^

.. automodule:: eSCAPE.flow.surfprocplex
    :members:

Pit filling
--------------------

Parallel priority-flood algorithm for unstructured mesh based on `Barnes, 2016 <https://arxiv.org/abs/1606.06204>`_

.. automodule:: eSCAPE.pit
    :members:

Depression algorithm
^^^^^^^^^^^^^^^^^^^^^

.. automodule:: eSCAPE.pit.pitfillunst
    :members:
