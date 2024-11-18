Installation
============

This section describes how to install `HTSinfer`, either from :ref:`Bioconda<Install from Bioconda>` (recommended for most users) or :ref:`GitHub<Install from GitHub>` (for experienced users and developers).

Both installation methods require the popular package and environment manager `Conda <https://docs.conda.io/en/latest/miniconda.html>`_ to be available on your system.

.. note::

   For improved performance, we recommend using Conda's drop-in replacement `Mamba <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_ instead. 
   In that case, simply replace :code:`conda` with :code:`mamba` in all of the commands below.

Install from Bioconda
---------------------

The easiest and quickest installation method is to install the latest `HTSinfer` release from `Bioconda <https://anaconda.org/bioconda/htsinfer>`_. 

To install `HTSinfer` in your currently active environment, run:

.. code-block:: bash

   conda install bioconda::htsinfer

.. note::

   We do not recommend installing HTSinfer into the default (:code:`base`) environment. 
   Instead, activate a project- or tool-specific environment with :code:`conda activate MY_CONDA_PROJECT`. 
   If you do not have such an environment, just create one with :code:`conda create --name MY_CONDA_PROJECT`. 
   :code:`MY_CONDA_PROJECT` can be any name you like (with some constraints imposed by Conda), e.g., :code:`htsinfer`.

Install from GitHub
-------------------

If you wish to use the latest unreleased version of the tool or if you would like to contribute to its further development,
please install `HTSinfer` from `GitHub <https://github.com/zavolanlab/htsinfer>`_.

First clone the repository:

.. code-block:: bash

   git clone https://github.com/zavolanlab/htsinfer.git
   cd htsinfer

Then, install `HTSinfer` and its dependencies via `Conda <https://docs.conda.io/en/latest/miniconda.html>`_:

.. code-block:: bash

   conda env create --file environment.yml

Alternatively, to include `HTSinfer`'s optional development dependencies, run the following instead:

.. code-block:: bash

   conda env create --file environment-dev.yml

After the installation is complete, activate the :code:`htsinfer` Conda environment with:

.. code-block:: bash

   conda activate htsinfer

If you have installed the development dependencies, you can verify that `HTSinfer` was installed correctly by executing the tests shipped with the package:

.. code-block:: bash

   python -m pytest
