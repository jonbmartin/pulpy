import math
import unittest

import numpy as np
import numpy.testing as npt

import pulpy.grad.optim as optim

if __name__ == "__main__":
    unittest.main()
    
class TestTrajGradOptim(unittest.TestCase):
    def test_min_gradient(self):
        t = np.linspace(0, 1, 1000)
        kx = np.sin(2.0 * math.pi * t)
        ky = np.cos(2.0 * math.pi * t)
        kz = t
        k = np.stack((kx, ky, kz), axis=-1)

        (g, k, s, t) = optim.min_time_gradient(
            k, 0.0, 0.0, gmax=4, smax=15, dt=4e-3, gamma=4.257
        )

        npt.assert_almost_equal(np.max(t), 0.916, decimal=4)
