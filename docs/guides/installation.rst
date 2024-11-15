Installation
============

This section describes how to install `HTSinfer`.

Install using Conda
-------------------

The easiest and quickest installation method is via `Mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_ or `Conda <https://docs.conda.io/en/latest/miniconda.html>`_.
`HTSinfer` is available as part of the `Bioconda <https://anaconda.org/bioconda/htsinfer>`_ channel.

To create a new Conda environment with `HTSinfer` and its dependencies installed, run:

.. code-block:: bash

   mamba create --name htsinfer bioconda::htsinfer


Then, activate the `htsinfer` Conda environment with:

.. code-block:: bash

   conda activate htsinfer


To install `HTSinfer` in your current environment, run:

.. code-block:: bash

   mamba install bioconda::htsinfer

Install from GitHub
-------------------

If you would like to contribute to the development of `HTSinfer`, or wishing to use unreleased versions, you can install `HTSinfer` from the `GitHub <https://github.com/zavolanlab/htsinfer>`_.
First clone the repository and install the dependencies via `Conda <https://docs.conda.io/en/latest/miniconda.html>`_:

.. code-block:: bash

   git clone https://github.com/zavolanlab/htsinfer.git
   cd htsinfer
   mamba env create --file environment.yml
   # Alternatively, to install with development dependencies,
   # run the following instead
   mamba env create --file environment-dev.ymls

After the installation is complete, activate the :code:`htsinfer` Conda environment with:

.. code-block:: bash

   conda activate htsinfer

Verify the Installation (Optional)
----------------------------------

If you have installed the development or testing dependencies, you can verify that `HTSinfer` was installed correctly by executing the tests shipped with the package:

.. code-block:: bash

   python -m pytest
