Functions (`pulpy`)
===================
.. automodule::
   pulpy
		
Adiabatic Pulse Design Functions
--------------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.adiabatic.bir4
    pulpy.adiabatic.hypsec
    pulpy.adiabatic.wurst
    pulpy.adiabatic.goia_wurst
    pulpy.adiabatic.bloch_siegert_fm

B1-Selective Pulse Design Functions
-----------------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.b1sel.dz_b1_rf
    pulpy.b1sel.dz_b1_gslider_rf
    pulpy.b1sel.dz_b1_hadamard_rf

RF Linear Operators
--------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.linop.PtxSpatialExplicit

Pulse Multibanding Functions
----------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.multiband.mb_rf
    pulpy.multiband.dz_pins

Optimal Control Design Functions
--------------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.optcont.blochsim
    pulpy.optcont.deriv

Parallel Transmit Pulse Designers
---------------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.ptx.stspa
    pulpy.ptx.stspk

RF Shimming Functions
--------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.shim.calc_shims
    pulpy.shim.init_optimal_spectral
    pulpy.shim.init_circ_polar

RF Pulse Simulation
--------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.sim.abrm
    pulpy.sim.abrm_nd
    pulpy.sim.abrm_hp
    pulpy.sim.abrm_ptx

SLR Pulse Design Functions
--------------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.slr.dzrf
    pulpy.slr.root_flip
    pulpy.slr.dz_gslider_rf
    pulpy.slr.dz_gslider_b
    pulpy.slr.dz_hadamard_b
    pulpy.slr.dz_recursive_rf

Trajectory and Gradient Design Functions
----------------------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.trajgrad.min_trap_grad
    pulpy.trajgrad.trap_grad
    pulpy.trajgrad.spiral_varden
    pulpy.trajgrad.spiral_arch
    pulpy.trajgrad.epi
    pulpy.trajgrad.rosette
    pulpy.trajgrad.stack_of
    pulpy.trajgrad.spokes_grad
    pulpy.trajgrad.traj_complex_to_array
    pulpy.trajgrad.traj_array_to_complex
    pulpy.trajgrad.min_time_gradient

RF Utility
--------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.util.dinf

I/O
--------------------------
.. autosummary::
    :toctree: generated
    :nosignatures:

    pulpy.io.siemens_rf
    pulpy.io.signa
    pulpy.io.ge_rf_params
    pulpy.io.philips_rf_params
