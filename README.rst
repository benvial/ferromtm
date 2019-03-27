ferromtm
==============================

.. inclusion-marker-do-not-remove

Coupled model and homogenization of ferroelectric-dielectric metamaterials.
---------------------------------------------------------------------------

This repository provides the codes to run and postprocess the data for the
results obtained in this research project.

Requirements
++++++++++++

- Python 3
- Gmsh_
- GetDP_
- make

Installation
++++++++++++

First clone this repository:

.. code-block:: bash

  git clone https://github.com/benvial/ferromtm.git

Then create, activate the environment and test it:

.. code-block:: bash

    cd ferromtm
    make env
    source activate ferromtm
    make testenv



Finally install the required packages:

.. code-block:: bash

  make req


Alternatively, one can use this docker recipe_:

.. code-block:: bash

   docker run -t --rm benvial/ferromtm:latest

.. _Gmsh: http://www.gmsh.info/
.. _GetDP: http://www.getdp.info/
.. _recipe: https://hub.docker.com/r/benvial/ferromtm
