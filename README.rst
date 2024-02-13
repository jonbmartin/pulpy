PulPy: Pulses in Python
=======================

.. image:: ../docs/figures/pulpy_logo_v2.png
  :width: 250
  :alt: Alternative text


`Source Code <https://github.com/jonbmartin/pulpy>`_ | `Documentation <https://pulpy.readthedocs.io>`_ |

.. image:: https://readthedocs.org/projects/pulpy/badge/?version=latest
    :target: https://pulpy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Description
-----------
Welcome to PulPy! PulPy is a python package for RF pulse and gradient pulse design for MRI.

PulPy is the successor package to `SigPy.RF <https://github.com/jonbmartin/sigpy-rf>`_, which is a sub package for RF pulse
design nested inside of the `SigPy <https://github.com/mikgroup/sigpy>`_ Python package for signal processing and image reconstruction.
PulPy is designed to replicate the functionality of SigPy.RF, but is able to stand alone as its own Python package for 
those with less interest in the other features of SigPy, such as iterative image reconstruction. It also places greater emphasis on 
MRI gradient waveform design, which is necessary for many classes of RF pulses but also has broader applicability.

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