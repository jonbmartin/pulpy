PulPy: Pulses in Python
=======================

.. image:: ../docs/figures/pulpy_logo_v2.png
  :width: 250
  :alt: Alternative text


`Source Code <https://github.com/jonbmartin/pulpy>`_ | `Documentation <https://pulpy.readthedocs.io>`_ |

.. image:: https://readthedocs.org/projects/pulpy/badge/?version=latest
    :target: https://pulpy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
    
## Description
TODO: this is a placeholder 

A brief description of what this project does and who it's for.

Installation
------------

PulPy requires Python version >= 3.8. The core module depends on ``numba``, ``numpy``, ``PyWavelets``, ``scipy``, and ``tqdm``.

Via ``pip``
***********

PulPy can be installed through ``pip``::
	
    pip install pulpy

Developer Installation
***************************

If you want to contribute to the PulPy source code, we recommend you install it with ``pip`` in editable mode::

	cd /path/to/pulpy
	pip install -e .
	
To run tests and contribute, we recommend installing the following packages::

	pip install coverage ruff sphinx sphinx_rtd_theme black isort

and run the script ``run_tests.sh``.