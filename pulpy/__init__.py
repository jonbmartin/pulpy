"""The core module contains functions and classes for RF pulse design.

It contains functions to design a variety of RF pulses for MRI, such as SLR,
adiabatic, parallel transmit, multibanded, and others. The submodule also
includes other functions to assist with pulse design, such as I/O functions,
trajectory/gradient designers, and Bloch simulators.

"""
from pulpy import (adiabatic, b1sel, io, multiband, optcont, ptx, shim, sim,
                   slr, trajgrad, util)
from pulpy.adiabatic import *  # noqa
from pulpy.b1sel import *  # noqa
from pulpy.io import *  # noqa
from pulpy.linop import *  # noqa
from pulpy.multiband import *  # noqa
from pulpy.optcont import *  # noqa
from pulpy.ptx import *  # noqa
from pulpy.shim import *  # noqa
from pulpy.sim import *  # noqa
from pulpy.slr import *  # noqa
from pulpy.trajgrad import *  # noqa
from pulpy.util import *  # noqa

__all__ = ["linop"]
__all__.extend(adiabatic.__all__)
__all__.extend(b1sel.__all__)
__all__.extend(io.__all__)
__all__.extend(multiband.__all__)
__all__.extend(optcont.__all__)
__all__.extend(ptx.__all__)
__all__.extend(sim.__all__)
__all__.extend(shim.__all__)
__all__.extend(slr.__all__)
__all__.extend(trajgrad.__all__)
__all__.extend(util.__all__)
