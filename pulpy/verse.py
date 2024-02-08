# -*- coding: utf-8 -*-
"""RF and Gradient VERSE-related functions.

"""

import numpy as np
from scipy.interpolate import CubicSpline

__all__ = ["flatten_bsse"]


def flatten_bsse(rfin, dtin, dtout=None):
    """Convert a Bloch-Siegert Selective Excitation (BSSE) pulse to have flat-amplitude.

    Args:
        rfin (complex vector): input BSSE RF samples
        dtin (float): intput RF dwell time
        dtout (float): desired output RF dwell time

    Returns:
        **rfout** (*array*): complex flattend BSSE waveform

    """

    if dtout is None:
        dtout = dtin

    rfin = rfin.flatten()

    # first remove any zero points
    rfin = rfin[rfin != 0]

    # determine effective source time points,
    # assuming dt * rfin has a flat amplitude of 1
    dteff = np.abs(rfin) * dtin
    tin = np.cumsum(dteff)

    # get an interpolator
    cs = CubicSpline(tin, rfin / np.abs(rfin))

    # get the vector of output points
    Nout = np.ceil(tin[-1] / dtout)
    tout = np.arange(Nout) * dtout

    rfout = cs(tout)

    return rfout
