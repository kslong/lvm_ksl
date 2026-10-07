#!/usr/bin/env python
# coding: utf-8
'''
                    Space Telescope Science Institute

Synopsis:

    Fit and apply an empirical, wavelength-dependent RELATIVE line-spread-
    function (LSF) kernel between two LVM spectra -- e.g. a background
    ("sky") spectrum built from some fibers and the spectrum of another
    fiber (or a median of fibers) -- so that the sharper of the two can be
    broadened to match the other before subtracting the sky lines.

    ADAPTED FROM Ivan Katkov's sky_decomp.lsf_surface_iterative (lvmsky
    repository, skysub/sky_decomp/lsf_surface_iterative.py; reference
    commit e30b730), which fits a constrained 11-tap kernel whose taps are
    B-splines in wavelength.  This is an independent, simplified
    re-implementation that depends only on numpy/scipy; see Notes for what
    differs.

Description::

    Model, per spectrograph arm (b, r, z; boundaries ARMS):

        target(lambda_i) = sum_t k_t(lambda_i) * source(i - t)  +  P(lambda_i)

    t = -5..5 (N_TAPS pixel offsets), k_t(lambda) = sum_j c_tj B_j(lambda)
    with B_j cubic B-splines in wavelength (a single constant kernel, j = 1,
    in the b arm, which has few isolated sky lines -- fitted only redward of
    B_FIT_MIN), and P a low-order Legendre polynomial (a nuisance
    background).  The source is never shifted across an arm boundary.

    Constraints (the B-splines are non-negative and sum to one, so these
    make every evaluated kernel non-negative and normalized):
        c_tj >= 0;   sum_t c_tj = 1 for each j;
        optionally (unimodal=True, default) single-peaked at the central
        tap: c_tj non-decreasing up to t = 0 and non-increasing after.
    Regularization: a roughness penalty on the second difference of each
    tap's coefficients across j, and a weak pull toward a narrow Gaussian
    default kernel where the data do not constrain it.  Solved as a
    quadratic program with scipy.optimize.minimize (SLSQP).

    Because a kernel only broadens, two_sided_match() fits both directions
    (source -> target and target -> source) and at each wavelength
    broadens whichever spectrum is sharper (the direction whose kernel is
    wider), with the decision smoothed along wavelength.  Kernels conserve
    flux: throughput (flat-field) differences are NOT absorbed.

Primary routines::

    fit_kernel        fit one relative LSF kernel surface (source -> target)
    apply_kernel      convolve a spectrum with a fitted kernel surface
    kernel_moments    kernel centroid and width (pixels) at each wavelength
    two_sided_match   fit both directions and broaden the sharper spectrum

Notes:

    Differences from lsf_surface_iterative: no iterative refinement of a
    full sky decomposition (this only relates two observed spectra); knots
    are uniform or explicit (no information-quantile placement); a general
    SLSQP quadratic program instead of clarabel; simplified priors.
    Runs in the ksl environment.

History::

    260930 ksl Coding begun (Vela background-sky tests, test_sky/background),
        adapting Ivan Katkov's lsf_surface_iterative approach so that the
        LSF matching is independent of the lvmsky repository.  Moved to
        py_progs/ on branch sky_patch.

'''
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize
from scipy.ndimage import uniform_filter1d
from numpy.polynomial import legendre

N_TAPS = 11
TAPS = np.arange(-(N_TAPS // 2), N_TAPS // 2 + 1)
ARMS = (('B', None, 5787.0), ('R', 5787.0, 7454.0), ('Z', 7454.0, None))
B_FIT_MIN = 5500.0
DEFAULT_ARM_CONFIG = {
    'B': dict(n_basis=1, knots=None),
    'R': dict(n_basis=6, knots=None),
    'Z': dict(n_basis=12, knots=(7700.0, 7750.0, 7800.0, 7850.0, 7900.0, 8000.0, 8300.0, 9000.0)),
}
BG_DEGREE = 3
ROUGHNESS = 0.01          # roughness penalty, relative to the mean data curvature
PRIOR = 1.0e-4            # pull toward the default kernel, relative to the mean data curvature
DEFAULT_SIGMA_PIX = 0.3   # width of the narrow default kernel


def arm_mask(wave, lo, hi):
    return ((wave >= lo) if lo else np.ones(wave.size, bool)) & ((wave < hi) if hi else np.ones(wave.size, bool))


def _default_kernel():
    k = np.exp(-0.5 * (TAPS / DEFAULT_SIGMA_PIX) ** 2)
    return k / k.sum()


def _basis(x, n_basis, knots, lo, hi):
    '''B-spline design matrix (npix, n_basis) on [lo, hi]; constant if n_basis == 1'''
    if n_basis == 1:
        return np.ones((x.size, 1))
    k = 3
    interior = np.asarray(knots, float) if knots is not None else np.linspace(lo, hi, n_basis - k + 1)[1:-1]
    if interior.size != n_basis - k - 1:
        raise ValueError('need n_basis - 4 interior knots, got %d for n_basis %d' % (interior.size, n_basis))
    t = np.concatenate([[lo] * (k + 1), interior, [hi] * (k + 1)])
    return BSpline.design_matrix(np.clip(x, lo, hi - 1e-9 * (hi - lo)), t, k).toarray()


def _shifted(src, offset):
    '''src(i - offset), zero outside'''
    out = np.zeros_like(src)
    if offset > 0:
        out[offset:] = src[:-offset]
    elif offset < 0:
        out[:offset] = src[-offset:]
    else:
        out[:] = src
    return out


def fit_kernel(wave, target, source, ivar, arm_config=None, unimodal=True):
    '''
    Fit the relative kernel surface mapping source onto target (both are
    line spectra, i.e. continuum already removed; pixels to ignore have
    ivar = 0).

    Returns:
        dict with 'surface' (npix, N_TAPS) kernel at every pixel (rows sum
        to 1), 'coef' / 'range' per arm, and 'status' per arm
    '''
    arm_config = arm_config or DEFAULT_ARM_CONFIG
    wave = np.asarray(wave, float)
    surface = np.tile(_default_kernel(), (wave.size, 1))
    out = dict(surface=surface, coef={}, range={}, status={})
    for arm, lo, hi in ARMS:
        inarm = arm_mask(wave, lo, hi)
        cfg = arm_config[arm]
        src_arm = np.where(inarm, np.nan_to_num(source), 0.0)
        fit = inarm & (ivar > 0) & np.isfinite(target)
        if arm == 'B':
            fit &= wave >= B_FIT_MIN
        if fit.sum() < 50:
            out['status'][arm] = 'too few pixels; default kernel'
            continue
        w_lo, w_hi = wave[inarm].min(), wave[inarm].max()
        J = cfg['n_basis']
        Bfit = _basis(wave[fit], J, cfg['knots'], w_lo, w_hi)
        # put target and source on an O(1) scale with the SAME factor (the
        # kernel relates them, so both must be scaled alike), and scale the
        # weights to match; otherwise the unit-amplitude background columns
        # swamp the kernel columns and the solver stops at its start point
        scale = np.nanmax(np.abs(target[fit])) or 1.0
        src_s = src_arm / scale
        y_s = np.nan_to_num(target[fit]) / scale
        sw = np.sqrt(ivar[fit]) * scale
        sw = sw / np.median(sw[sw > 0])       # only relative weights matter
        # design: kernel columns (tap-major: index t*J + j), then background columns
        cols = [Bfit * _shifted(src_s, int(o))[fit, None] for o in TAPS]
        A_k = np.concatenate(cols, axis=1)
        xl = 2 * (wave[fit] - w_lo) / (w_hi - w_lo) - 1
        A_b = legendre.legvander(xl, BG_DEGREE)
        A = np.concatenate([A_k, A_b], axis=1)
        Aw = A * sw[:, None]
        yw = y_s * sw
        H = Aw.T @ Aw
        g = Aw.T @ yw
        nk = N_TAPS * J
        mean_curv = np.mean(np.diag(H)[:nk]) or 1.0
        # regularization on the kernel coefficients
        f0 = np.repeat(_default_kernel(), J)
        R = np.eye(nk) * PRIOR * mean_curv
        rhs = PRIOR * mean_curv * f0
        if J >= 3:
            D = np.diff(np.eye(J), 2, axis=0)
            DtD = D.T @ D
            for t in range(N_TAPS):
                R[t * J:(t + 1) * J, t * J:(t + 1) * J] += ROUGHNESS * mean_curv * DtD
        H[:nk, :nk] += R
        g[:nk] += rhs
        # normalize by the kernel block's own curvature
        H /= mean_curv
        g /= mean_curv

        def obj(x):
            return 0.5 * x @ H @ x - g @ x

        def jac(x):
            return H @ x - g

        cons = []
        # sum over taps = 1 for each basis function
        E = np.zeros((J, nk + A_b.shape[1]))
        for j in range(J):
            E[j, j:nk:J] = 1.0
        cons.append(dict(type='eq', fun=lambda x, E=E: E @ x - 1.0, jac=lambda x, E=E: E))
        if unimodal:
            rows = []
            c0 = N_TAPS // 2
            for j in range(J):
                for t in range(N_TAPS - 1):
                    r = np.zeros(nk + A_b.shape[1])
                    if t < c0:      # rising toward the centre: c[t+1] - c[t] >= 0
                        r[(t + 1) * J + j], r[t * J + j] = 1.0, -1.0
                    else:           # falling after it: c[t] - c[t+1] >= 0
                        r[t * J + j], r[(t + 1) * J + j] = 1.0, -1.0
                    rows.append(r)
            U = np.array(rows)
            cons.append(dict(type='ineq', fun=lambda x, U=U: U @ x, jac=lambda x, U=U: U))
        bounds = [(0.0, None)] * nk + [(None, None)] * A_b.shape[1]
        x0 = np.concatenate([f0, np.zeros(A_b.shape[1])])
        res = minimize(obj, x0, jac=jac, bounds=bounds, constraints=cons, method='SLSQP',
                       options=dict(maxiter=3000, ftol=1e-12))
        coef = np.clip(res.x[:nk].reshape(N_TAPS, J), 0, None)
        Ball = _basis(wave[inarm], J, cfg['knots'], w_lo, w_hi)
        surf = Ball @ coef.T
        surf /= surf.sum(axis=1, keepdims=True)
        surface[inarm] = surf
        out['coef'][arm] = coef
        out['range'][arm] = (w_lo, w_hi)
        out['status'][arm] = 'ok' if res.success else 'SLSQP: ' + res.message
    return out


def apply_kernel(wave, spec, surface):
    '''convolve spec with a kernel surface, never mixing pixels across arm boundaries'''
    spec = np.nan_to_num(np.asarray(spec, float))
    out = np.zeros_like(spec)
    for arm, lo, hi in ARMS:
        inarm = arm_mask(wave, lo, hi)
        s_arm = np.where(inarm, spec, 0.0)
        acc = np.zeros_like(spec)
        for k, o in enumerate(TAPS):
            acc += surface[:, k] * _shifted(s_arm, int(o))
        out[inarm] = acc[inarm]
    return out


def kernel_moments(surface):
    '''centroid and standard deviation (pixels) of the kernel at each pixel'''
    c = surface @ TAPS
    s = np.sqrt(np.maximum(np.sum(surface * (TAPS[None, :] - c[:, None]) ** 2, axis=1), 0.0))
    return c, s


def two_sided_match(wave, target, source, ivar, fit_source=None, arm_config=None, unimodal=True, smooth_pix=51):
    '''
    Fit kernels in both directions and, at each wavelength, broaden
    whichever of (target, source) is sharper.

    Parameters:
        target, source: line spectra (continuum removed)
        ivar: weights; 0 for pixels to ignore (nebular lines, bad regions)
        fit_source: optional version of source used only for FITTING the
            kernels (e.g. scaled per arm so a throughput difference does
            not distort the kernel shape); the kernels are applied to source
        smooth_pix: smoothing length for the which-is-sharper decision

    Returns:
        dict with target_matched, source_matched, source_sharper (bool per
        pixel), k1 (source -> target fit), k2 (target -> source fit),
        sigma1, sigma2 (kernel widths, pixels)
    '''
    fs = source if fit_source is None else fit_source
    k1 = fit_kernel(wave, target, fs, ivar, arm_config, unimodal)
    k2 = fit_kernel(wave, fs, target, ivar, arm_config, unimodal)
    _, s1 = kernel_moments(k1['surface'])
    _, s2 = kernel_moments(k2['surface'])
    src_sharper = uniform_filter1d(s1 - s2, size=smooth_pix) >= 0
    tgt_m = np.where(src_sharper, target, apply_kernel(wave, target, k2['surface']))
    src_m = np.where(src_sharper, apply_kernel(wave, source, k1['surface']), source)
    return dict(target_matched=tgt_m, source_matched=src_m, source_sharper=src_sharper,
                k1=k1, k2=k2, sigma1=s1, sigma2=s2)
