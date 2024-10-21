Installation
============

This section describes how to install HTSinfer using Conda.

Clone the Repository and Install Dependencies
---------------------------------------------

To install HTSinfer, first clone the repository and install the dependencies via `Conda <https://docs.conda.io/en/latest/miniconda.html>`_:

.. code-block:: bash

   git clone https://github.com/zavolanlab/htsinfer.git
   cd htsinfer
   conda env create --file environment.yml
   # Alternatively, to install with development dependencies,
   # run the following instead
   conda env create --file environment-dev.ymls

.. note::

   Creating the environment may take some time. It is strongly recommended to install `Mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_ and replace ``conda`` with ``mamba`` in the previous commands for faster installation.

Activate the Conda Environment
------------------------------

After the installation is complete, activate the `htsinfer` Conda environment with:

.. code-block:: bash

   conda activate htsinfer

Verify the Installation (Optional)
----------------------------------

If you have installed the development or testing dependencies, you can verify that HTSinfer was installed correctly by executing the tests shipped with the package:

.. code-block:: bash

   python -m pytest
