PulPy: Pulses in Python
=====

.. image:: docs/figures/pulpy_logo_v2.png
  :width: 250
  :alt: Alternative text


`Source Code <https://github.com/mikgroup/sigpy>`_ | `Documentation <https://sigpy.readthedocs.io>`_ |


## Description
TODO: this is a placeholder 

A brief description of what this project does and who it's for.

Installation
------------

SigPy requires Python version >= 3.5. The core module depends on ``numba``, ``numpy``, ``PyWavelets``, ``scipy``, and ``tqdm``.

Via ``pip``
***********

SigPy can be installed through ``pip``::
	
    pip install sigpy

Developer Installation
***************************

If you want to contribute to the PulPy source code, we recommend you install it with ``pip`` in editable mode::

	cd /path/to/pulpy
	pip install -e .
	
To run tests and contribute, we recommend installing the following packages::

	pip install coverage ruff sphinx sphinx_rtd_theme black isort

and run the script ``run_tests.sh``.