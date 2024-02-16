PulPy: Pulses in Python
=======================

.. image:: ../docs/figures/pulpy_logo_v2.png
  :width: 250
  :alt: Alternative text


`Source Code <https://github.com/jonbmartin/pulpy>`_ | `Documentation <https://pulpy.readthedocs.io>`_

.. image:: https://readthedocs.org/projects/pulpy/badge/?version=latest
    :target: https://pulpy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Description
-----------
Welcome to PulPy! PulPy is a python package for RF pulse and gradient pulse design for MRI.

PulPy is the successor package to `SigPy.RF <https://github.com/jonbmartin/sigpy-rf>`_, a sub package for RF pulse
design nested inside of the `SigPy <https://github.com/mikgroup/sigpy>`_ package for signal processing and image reconstruction.
PulPy is designed to replicate  functionality of SigPy.RF, but is able to stand alone as an independent Python package. It also places greater emphasis on 
MRI gradient waveform design, which is necessary for many classes of RF pulses but also has broader applicability. PulPy current relies on SigPy as a 
dependency and takes advantage of its' abstraction classes. This dependency will be reduced in future releases. 

Installation
------------

PulPy requires Python version >= 3.8. The core module depends on ``numba``, ``numpy``, ``sigpy``, ``PyWavelets``, ``scipy``, and ``tqdm``.

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

If you make modifications and would like to update the version number across the release, run bumpversion 
(where argument size is one of 'major', 'minor', 'patch')::

  bumpversion --current-version 1.8.1 size

Contact and Contribution
------------------------
We welcome feedback on this project! It is a work in project, so please report bugs and issues on 
GitHub. We also encourage you to contribute additional pulse design tools. Point of contact: jonathan.bach.martin@vumc.org. 