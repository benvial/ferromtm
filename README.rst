ferromtm
==============================

.. inclusion-marker-do-not-remove

Coupled model and homogenization of ferroelectric-dielectric metamaterials.
---------------------------------------------------------------------------

This repository provides the codes to run, postprocess the results and compile the
article related to this research project.

Requirements
++++++++++++

- Python 3
- Gmsh_
- GetDP_
- make
- pdflatex and bibtex (for the paper)
- xelatex and biblatex/biber (for the poster)

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


Running the code
++++++++++++++++

Make the data (this will take a while, if you have access to a cluster with
Sun Grid Engine this should be parallelized automatically)

.. code-block:: bash

  make results


Postprocessing and plots
++++++++++++++++++++++++

This will postprocess some of the data

.. code-block:: bash

  make postpro


Plotting
++++++++++++++++++++++++

Generate plots with:

.. code-block:: bash

  make plots



Article
+++++++

Run pdflatex and generate the pdf paper

.. code-block:: bash

  make paper


Poster
+++++++

Run xelatex and generate the pdf poster

.. code-block:: bash

  make poster
