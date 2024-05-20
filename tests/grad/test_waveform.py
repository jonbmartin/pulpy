import math
import unittest

import numpy as np
import numpy.testing as npt

import pulpy.grad.waveform as waveform

if __name__ == "__main__":
    unittest.main()


class TestTrajGrad(unittest.TestCase):
    def test_trap_grad(self):
        dt = 4e-6  # s
        area = 200 * dt
        dgdt = 18000  # g/cm/s
        gmax = 2  # g/cm

        trap, _ = waveform.trap_grad(area, gmax, dgdt, dt)

        npt.assert_almost_equal(area, np.sum(trap) * dt, decimal=3)
        npt.assert_almost_equal(gmax, np.max(trap), decimal=1)

    def test_min_trap_grad(self):
        dt = 4e-6  # s
        area = 200 * dt
        dgdt = 18000  # g/cm/s
        gmax = 2  # g/cm

        trap, _ = waveform.min_trap_grad(area, gmax, dgdt, dt)

        npt.assert_almost_equal(area, np.sum(trap) * dt, decimal=3)
        npt.assert_almost_equal(gmax, np.max(trap), decimal=1)
