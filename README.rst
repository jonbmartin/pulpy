PulPy: Pulses in Python
=======================

.. image:: /docs/figures/pulpy_logo_v2.png
   :align: center
   :width: 250
   :alt: PulPy logo


`Source Code <https://github.com/jonbmartin/pulpy>`_ | `Documentation <https://pulpy.readthedocs.io>`_ | `Demo Code <https://github.com/jonbmartin/pulpy-tutorials>`_

.. image:: https://readthedocs.org/projects/pulpy/badge/?version=latest
    :target: https://pulpy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Description
-----------
Welcome to PulPy! PulPy is a Python package for RF pulse and gradient pulse design for MRI. PulPy focuses on the design of individual pulses rather than full pulse sequences (for this, see other packages, e.g. `PyPulSeq <https://github.com/imr-framework/pypulseq>`_).

PulPy is the successor package to `SigPy.RF <https://github.com/jonbmartin/sigpy-rf>`_, a sub package for RF pulse
design nested inside of the `SigPy <https://github.com/mikgroup/sigpy>`_ package for signal processing and image reconstruction.
PulPy is designed to replicate  functionality of SigPy.RF, but is able to stand alone as an independent Python package. It also places greater emphasis on 
MRI gradient waveform design, which is necessary for many classes of RF pulses but also has broader applicability. PulPy current relies on SigPy as a 
dependency and takes advantage of its' abstraction classes. This dependency will be reduced in future releases. 

Installation
------------

PulPy requires Python version >= 3.9 and has been tested through 3.12. The core module depends on ``numba``, ``numpy``, ``sigpy``, ``PyWavelets``, ``scipy``, and ``tqdm``.

Via ``pip``
***********

PulPy can be installed through ``pip``::
	
    pip install pulpy

Developer Installation
***************************

If you want to contribute to the PulPy source code, we recommend you install it with ``pip`` in editable mode::

	cd /path/to/pulpy
	pip install -e .
	
Tests are currently run through GitHub Actions on submission of pull requests or on commits to the main branch. 


Getting Started
**********************
To begin using `pulpy`, import the package in your Python script. For demo purposes, we'll also import the SigPy package's plotting functions:

.. code-block:: python

   import pulpy as pp    	# import full package
   import sigpy.plot as pl      # import a plotting function

1) RF Pulse Design and Simulation
**************************************
`pulpy` allows you to easily design and simulate RF pulses. Here's an example of an SLR pulse design with multibanding: 

1a) define RF pulse parameters 

.. code-block:: python
	
	tb = 8 		# RF pulse time-bandwidth product
	N = 512 	# number of timepoints to design
	d1 = 0.01 	# magnetization passband ripple level
	d2 = 0.01 	# magnetization stopband ripple level
	p_type = 'ex'   # RF pulse type - a 90 degree excitation pulse
	f_type = 'ls'   # filter type for SLR design - using a least squares filter

1b) design and plot the RF pulse

.. code-block:: python

	pulse = pp.rf.slr.dzrf(N, tb, p_type, f_type, d1, d2, True)
	pl.LinePlot(pulse,mode='r')     # plot the real component of the RF pulse

.. image:: /docs/figures/slr_pulse.png
   :align: center
   :width: 300

1c) multiband the single-band RF pulse to excite multiple slices simultaneously

.. code-block:: python

	n_bands = 3              # design to excite 3 bands of magnetizaztion
	phs_type = 'phs_mod'     # 'phsMod', 'ampMod', or 'quadMod' - the method of designing the pulse phases
	band_sep = 5*tb          # separate by 5 slice widths
	mb_pulse = pp.rf.multiband.mb_rf(pulse, n_bands, band_sep, phs_type)
	pl.LinePlot(mb_pulse)

.. image:: /docs/figures/multiband_pulse.png
   :align: center
   :width: 300

1d) simulate the transverse magnetization profile of both pulses. We do this by first calculating the Cayley-Klein parameters representing the rotation of the magnetization vector produced by the RF pulse (variables 'a' and 'b'). We then use the relationships in Pauly et. al. to convert this to the resulting excitation magnetization. 

.. code-block:: python

	[a, b] = pp.sim.abrm(pulse, np.arange(-20*tb, 20*tb, 40*tb/2000), True)
	Mxy_single_band = 2*np.multiply(np.conj(a), b)  # from Pauly et. al. IEEE TMI (1991). 
	[a, b] = pp.sim.abrm(mb_pulse, np.arange(-20*tb, 20*tb, 40*tb/2000), True)
	Mxy_multi_band = 2*np.multiply(np.conj(a), b)  # from Pauly et. al. IEEE TMI (1991). 
	pl.LinePlot(Mxy_single_band, title='single band excitation')
	pl.LinePlot(Mxy_multi_band, title='multi-band excitation')

.. image:: /docs/figures/single_band_excitation.png
   :align: center
   :width: 300
.. image:: /docs/figures/multiband_excitation.png
   :align: center
   :width: 300

1e) Export the RF pulse to GE format for use in a scanner. We will compute the important parameters
then write to .i file:

.. code-block:: python

	pp.ge_rf_params(pulse, dt=4e-6)   # prints out the most important GE parameters
	pp.signa(pulse, 'slr_ex')         # writes to .i file


2) Gradient Waveform Design and Optimization
************************************************
`pulpy` also has a variety of tools for designing gradient pulses. This ranges from simple trapezoids, the building block of many pulse sequences: 

.. code-block:: python

        dt = 4e-6  # s
        area = 200 * dt
        dgdt = 18000  # g/cm/s
        gmax = 2  # g/cm

        trap, _ = pp.grad.waveform.trap_grad(area, gmax, dgdt, dt)
        
        pl.LinePlot(trap, title='trapezoidal gradient')

.. image:: /docs/figures/trap_grad.png
   :align: center
   :width: 300

to more complex time-varying waveforms (e.g. spiral gradient waveform):

.. code-block:: python

        fov = 0.55    # imaging field of view [m]
        gts = 6.4e-6  # hardware dwell time [s]
        gslew = 190   # max. slew rate [mT/m/ms]
        gamp = 40     # max. amplitude [mT/m]
        R = 1         # degree of undersampling
        dx = 0.025    # resolution
        
        # construct a trajectory
        g, k, t, s = pp.grad.waveform.spiral_arch(fov / R, dx, gts, gslew, gamp)
        
        pl.LinePlot(np.transpose(g),mode='r', title='spiral gradient (1 axis plotted)')

.. image:: /docs/figures/spiral_waveform.png
   :align: center
   :width: 300

to a few tools for more advanced design (e.g. min-time-gradient designers, which modifies an existing trajectory to be
time-efficient): 

.. code-block:: python

	import math        
    
	t = np.linspace(0, 1, 1000)
	kx = np.sin(2.0 * math.pi * t)
	ky = np.cos(2.0 * math.pi * t)
	kz = t
	k = np.stack((kx, ky, kz), axis=-1)
	
	(g, k, s, t) = pp.grad.optim.min_time_gradient(
	    k, 0.0, 0.0, gmax=4, smax=15, dt=4e-3, gamma=4.257
	)

Documentation
**************
Documentation for PulPy is available at `ReadTheDocs <https://pulpy.readthedocs.io>`_.

A series of Jupyter notebooks have been developed that provide tutorials of several classes 
of pulse design at `the demo code repository <https://github.com/jonbmartin/pulpy-tutorials>`_.
Simply clone this repository, install Pulpy (and Jupyter notebook), and get started designing pulses! 

Contact and Contribution
*************************

We welcome feedback on this project! It is a work in project, so please report bugs and issues on 
GitHub. We also encourage you to contribute additional pulse design tools. Point of contact: jonathan.bach.martin@vumc.org. 
