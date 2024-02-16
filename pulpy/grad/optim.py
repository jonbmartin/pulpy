# -*- coding: utf-8 -*-
"""Utilities for optimizing existing gradient waveforms, such as M. Lustig's min-time-gradient designer
"""

import math

import numba as nb
import numpy as np
from scipy import integrate, interpolate

__all__ = [
    "min_time_gradient",
]


@nb.jit(nopython=True, cache=True)  # pragma: no cover
def runge_kutta(ds: float, st: float, kvals: np.ndarray, smax=None, gamma=4.257):
    r"""Runge-Kutta 4 for curve constrained

    Args:
        ds (float): spacing in arc length space
        st (float): output shape.
        kvals (array): 3 points of curve.
        smax (float): maximum slew
        gamma (float): gyromagnetic ratio

    Returns:
        float or None: step size dsdt or None
    """
    temp = gamma**2 * smax**2 - abs(kvals[0]) ** 2 * st**4
    if temp < 0.0:
        return None
    k1 = ds / st * math.sqrt(temp)

    temp = gamma**2 * smax**2 - abs(kvals[1]) ** 2 * (st + ds * k1 / 2) ** 4
    if temp < 0.0:
        return None
    k2 = ds / (st + ds * k1 / 2) * math.sqrt(temp)

    temp = gamma**2 * smax**2 - abs(kvals[1]) ** 2 * (st + ds * k2 / 2) ** 4
    if temp < 0.0:
        return None
    k3 = ds / (st + ds * k2 / 2) * math.sqrt(temp)

    temp = gamma**2 * smax**2 - abs(kvals[2]) ** 2 * (st + ds * k3) ** 4
    if temp < 0.0:
        return None
    k4 = ds / (st + ds * k3) * math.sqrt(temp)

    return k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6


#  Arc length code translated from matlab
#    (c) Michael Lustig 2005
#    modified 2006 and 2007
#    Rewritten in Python in 2020 by Kevin Johnson
def min_time_gradient(
    c: np.ndarray, g0=0, gfin=None, gmax=4, smax=15, dt=4e-3, gamma=4.257
):
    r"""
    Given a k-space trajectory c(n), gradient and slew constraints. This
    function will return a new parametrization that will meet these
    constraint while getting from one point to the other in minimum time.

    Args:
        c (array): Curve in k-space given in any parametrization [1/cm]
                        Nx3 real array
        g0 (float): Initial gradient amplitude (leave empty for g0 = 0)
        gfin (float): Gradient value at the end of the trajectory. If not
                        possible, the result would be the largest possible
                        ampltude. (Leave empty if you don't care to get
                        maximum gradient.)
        gmax (float): Maximum gradient [G/cm] (3.9 default)
        smax (float): Maximum slew [G/Cm/ms]  (14.5 default)
        dt (float): Sampling time interval [ms] (4e-3 default)
        gamma (float): Gyromagnetic ratio

    Returns:
        tuple: (g, k, s, t) tuple containing

        - **g** - (array): gradient waveform [G/cm]
        - **k** - (array): exact k-space corresponding to gradient g.
        - **s** - (array): slew rate [G/cm/ms]
        - **time** - (array):  sampled time

    References:
        Lustig M, Kim SJ, Pauly JM. A fast method for designing time-optimal
        gradient waveforms for arbitrary k-space trajectories. IEEE Trans Med
        Imaging. 2008;27(6):866-873. doi:10.1109/TMI.2008.922699
    """

    def sdotmax(cs: interpolate.CubicSpline, s: np.ndarray, gmax, smax, gamma=4.257):
        # [sdot, k, ] = sdotMax(PP, p_of_s, s, gmax, smax)
        #
        # Given a k-space curve C (in [1/cm] units), maximum gradient amplitude
        # (in G/cm) and maximum slew-rate (in G/(cm*ms)).
        # This function calculates the upper bound for the time parametrization
        # sdot (which is a non scaled max gradient constaint) as a function
        # of s.
        #
        #   cs      --  spline polynomial
        #   p_of_s  --  parametrization vs arclength
        #   s       --  arclength parametrization (0->1)
        #   gmax    --  maximum gradient (G/cm)
        #   smax    --  maximum slew rate (G/ cm*ms)
        #
        #   returns the maximum sdot (1st derivative of s) as a function of
        #   arclength s
        #   Also, returns curvature as a function of s and length of curve (L)
        #
        #  (c) Michael Lustig 2005
        #  last modified 2006

        # Absolute value of 2nd derivative in curve space using cubic splines
        cs2 = cs.derivative(2)  # spline derivative
        cs2_highres = cs2(s)  # evaluated along arc length
        k = np.linalg.norm(cs2_highres, axis=1)  # magnitude

        # calc I constraint curve (maximum gradient)
        sdot1 = gamma * gmax * np.ones_like(s)

        # calc II constraint curve (curve curvature dependent)
        sdot2 = np.sqrt(gamma * smax / (k + np.finfo(float).eps))

        # calc total constraint
        sdot = np.minimum(sdot1, sdot2)

        return sdot, k

    # Curve in arbitrary paramater space, cubic spline
    num_p = c.shape[0]
    p = np.linspace(0, 1, num_p, endpoint=True)
    cp = interpolate.CubicSpline(p, c, axis=0)

    # Integrate absolute value to find length and s(arc) vs p(paramater)
    cp1_spline = cp.derivative()
    p_highres = np.linspace(0, 1, num_p * 10)
    cp1_highres = cp1_spline(p_highres)
    ds_p = np.linalg.norm(cp1_highres, axis=1)

    # s vs p to enable conversion
    s_of_p = integrate.cumtrapz(ds_p, p_highres, initial=0)
    curve_length = s_of_p[-1]

    # decide ds and compute st for the first point
    stt0 = gamma * smax  # always assumes first point is max slew
    st0 = stt0 * dt / 8  # start at 1/8 the gradient for accuracy close to g=0
    s0 = st0 * dt
    ds = s0 / 4.0  # smaller step size for numerical accuracy
    ns = int(curve_length / ds)

    if g0 is None:
        g0 = 0

    # s is arc length at high resolution
    s = np.linspace(0, curve_length, ns, endpoint=True)

    # Cubic spline at s positions (s of p)
    cp_highres = cp(p_highres)
    cs = interpolate.CubicSpline(s_of_p, cp_highres, axis=0)

    # compute constraints (forbidden line curve)
    phi, k = sdotmax(cs, s, gmax, smax)

    # extend for the Runge-Kutte method
    k = np.pad(k, (0, 3), "constant", constant_values=(0,))

    # Get the start
    sta = np.zeros_like(s)
    sta[0] = min(g0 * gamma + st0, gamma * gmax)

    # solve ODE forward
    for n in range(1, s.shape[0]):
        kpos = n
        dstds = runge_kutta(ds, sta[n - 1], k[kpos : kpos + 4], smax)

        if dstds is None:
            sta[n] = phi[n]
        else:
            tmpst = sta[n - 1] + dstds
            sta[n] = min(tmpst, phi[n])

    stb = 0 * s
    if gfin is None:
        stb[-1] = sta[-1]
    else:
        stb[-1] = min(max(gfin * gamma, st0), gamma * gmax)

    # solve ODE backwards
    for n in range(s.shape[0] - 2, 0, -1):
        kpos_end = n  # to 0
        kpos = kpos_end + 3
        dstds = runge_kutta(ds, stb[n + 1], k[kpos : (kpos - 3) : -1], smax)

        if dstds is None:
            stb[n] = phi[n - 1]
        else:
            tmpst = stb[n + 1] + dstds
            stb[n] = min(tmpst, phi[n - 1])

    # Fix last point which is indexed a bit off
    n = 0
    kpos_end = n
    kpos = kpos_end + 3
    dstds = runge_kutta(ds, stb[n + 1], k[kpos::-1], smax)
    if dstds is None:
        stb[n] = phi[n * 2 - 1]
    else:
        tmpst = stb[n + 1] + dstds
        stb[n] = min(tmpst, phi[n - 1])

    # take the minimum of the curves
    ds = s[1] - s[0]
    st_of_s = np.minimum(sta, stb)

    # compute time
    t_of_s = integrate.cumtrapz(1.0 / st_of_s, initial=0) * ds

    t = np.arange(0, t_of_s[-1] + np.finfo(float).eps, dt)

    t_of_s = interpolate.CubicSpline(t_of_s, s)
    s_of_t = t_of_s(t)
    c = np.squeeze(cs(s_of_t))

    g = np.diff(c, axis=0, append=np.zeros((1, 3))) / gamma / dt
    g[-1, :] = g[-2, :] + g[-2, :] - g[-3, :]

    k = integrate.cumtrapz(g, t, initial=0, axis=0) * gamma

    s = np.diff(g, axis=0) / dt

    return g, k, s, t
