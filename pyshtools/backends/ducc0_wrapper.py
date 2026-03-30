"""
    Distinctly Useful Code Collection (DUCC)

    DUCC wrapper functions for use in pyshtools.
"""
import numpy as _np

try:
    import ducc0

    major, minor, patch = ducc0.__version__.split(".")
    if int(major) < 1 and int(minor) < 15:
        raise RuntimeError
except Exception:
    ducc0 = None

# setup a few required variables
if ducc0 is not None:
    import os as _os

    try:
        nthreads = int(_os.environ["OMP_NUM_THREADS"])
    except Exception:
        nthreads = 0


def _fixdtype(arr):
    return arr.astype(_np.float64, copy=False)


def set_nthreads(ntnew):
    global nthreads
    nthreads = ntnew


def available():
    return ducc0 is not None


def _nalm(lmax, mmax):
    return ((mmax + 1) * (mmax + 2)) // 2 + (mmax + 1) * (lmax - mmax)


# ducc0's only accepted conventions are
#   normalization = 'ortho'
#   csphase = 1
# so we need to make the required adjustments


def _get_norm(lmax, norm):
    if norm == 1:
        return _np.full(lmax + 1, _np.sqrt(4 * _np.pi))
    if norm == 2:
        return _np.sqrt(4 * _np.pi / (2 * _np.arange(lmax + 1) + 1.0))
    if norm == 3:
        return _np.sqrt(2 * _np.pi / (2 * _np.arange(lmax + 1) + 1.0))
    if norm == 4:
        return _np.ones(lmax + 1)
    raise RuntimeError("unsupported normalization")


def _rcilm2alm(cilm, lmax):
    alm = _np.empty((_nalm(lmax, lmax),), dtype=_np.complex128)
    alm[0: lmax + 1] = cilm[0, :, 0]
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        alm[ofs: ofs + lmax + 1 - m].real = cilm[0, m:, m]
        alm[ofs: ofs + lmax + 1 - m].imag = cilm[1, m:, m]
        ofs += lmax + 1 - m
    return alm


def _ralm2cilm(alm, lmax):
    cilm = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.float64)
    cilm[0, :, 0] = alm[0: lmax + 1].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        cilm[0, m:, m] = alm[ofs: ofs + lmax + 1 - m].real
        cilm[1, m:, m] = alm[ofs: ofs + lmax + 1 - m].imag
        ofs += lmax + 1 - m
    return cilm


def _apply_norm(alm, lmax, norm, csphase, reverse):
    lnorm = _get_norm(lmax, norm)
    if reverse:
        lnorm = 1.0 / lnorm
    alm[0: lmax + 1] *= lnorm[0: lmax + 1]
    lnorm *= _np.sqrt(2.0) if reverse else (1.0 / _np.sqrt(2.0))
    mlnorm = -lnorm
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        if csphase == 1:
            if m & 1:
                alm[ofs: ofs + lmax + 1 - m].real *= mlnorm[m:]
                alm[ofs: ofs + lmax + 1 - m].imag *= lnorm[m:]
            else:
                alm[ofs: ofs + lmax + 1 - m].real *= lnorm[m:]
                alm[ofs: ofs + lmax + 1 - m].imag *= mlnorm[m:]
        else:
            alm[ofs: ofs + lmax + 1 - m].real *= lnorm[m:]
            alm[ofs: ofs + lmax + 1 - m].imag *= mlnorm[m:]
        ofs += lmax + 1 - m
    if norm == 3:  # special treatment for unnormalized a_lm
        r = _np.arange(lmax + 1)
        fct = _np.ones(lmax + 1)
        ofs = lmax + 1
        if reverse:
            alm[0: lmax + 1] /= _np.sqrt(2)
            for m in range(1, lmax + 1):
                fct[m:] *= _np.sqrt((r[m:] + m) * (r[m:] - m + 1))
                alm[ofs: ofs + lmax + 1 - m] /= fct[m:]
                ofs += lmax + 1 - m
        else:
            alm[0: lmax + 1] *= _np.sqrt(2)
            for m in range(1, lmax + 1):
                fct[m:] *= _np.sqrt((r[m:] + m) * (r[m:] - m + 1))
                alm[ofs: ofs + lmax + 1 - m] *= fct[m:]
                ofs += lmax + 1 - m
    return alm


def _make_alm(cilm, lmax, norm, csphase):
    alm = _rcilm2alm(cilm, lmax)
    return _apply_norm(alm, lmax, norm, csphase, False)


def _extract_alm(alm, lmax, norm, csphase):
    _apply_norm(alm, lmax, norm, csphase, True)
    return _ralm2cilm(alm, lmax)


# ---- DH geometry ----

def _synthesize_DH(alm, lmax, extend, out):
    ducc0.sht.experimental.synthesis_2d(
        alm=alm.reshape((1, -1)),
        map=out[:, : out.shape[1] - extend].reshape(
            (1, out.shape[0], out.shape[1] - extend)
        ),
        spin=0,
        lmax=lmax,
        geometry="DH",
        nthreads=nthreads,
    )
    if extend:
        out[:, -1] = out[:, 0]
    return out


def _synthesize_DH_deriv1(alm, lmax, extend, out):
    ducc0.sht.experimental.synthesis_2d_deriv1(
        alm=alm.reshape((1, -1)),
        map=out[:, :, : out.shape[2] - extend],
        lmax=lmax,
        geometry="DH",
        nthreads=nthreads,
    )
    out[:, 0, :] = 0.0
    if extend:
        out[:, -1, :] = 0.0
        out[:, :, -1] = out[:, :, 0]
    return out


def _analyze_DH(map, lmax):
    alm = ducc0.sht.experimental.analysis_2d(
        map=map.reshape((1, map.shape[0], map.shape[1])),
        spin=0,
        lmax=lmax,
        geometry="DH",
        nthreads=nthreads,
    )
    return alm[0]


# ---- CC geometry ----

def _synthesize_CC(alm, lmax, extend, out):
    ducc0.sht.experimental.synthesis_2d(
        alm=alm.reshape((1, -1)),
        map=out[:, : out.shape[1] - extend].reshape(
            (1, out.shape[0], out.shape[1] - extend)
        ),
        spin=0,
        lmax=lmax,
        geometry="CC",
        nthreads=nthreads,
    )
    if extend:
        out[:, -1] = out[:, 0]
    return out


def _synthesize_CC_deriv1(alm, lmax, extend, out):
    ducc0.sht.experimental.synthesis_2d_deriv1(
        alm=alm.reshape((1, -1)),
        map=out[:, :, : out.shape[2] - extend],
        lmax=lmax,
        geometry="CC",
        nthreads=nthreads,
    )
    out[:, 0, :] = 0.0
    if extend:
        out[:, -1, :] = 0.0
        out[:, :, -1] = out[:, :, 0]
    return out


def _analyze_CC(map, lmax):
    alm = ducc0.sht.experimental.analysis_2d(
        map=map.reshape((1, map.shape[0], map.shape[1])),
        spin=0,
        lmax=lmax,
        geometry="CC",
        nthreads=nthreads,
    )
    return alm[0]


# ---- GLQ geometry ----

def _synthesize_GLQ(alm, lmax, extend, out):
    ducc0.sht.experimental.synthesis_2d(
        alm=alm.reshape((1, -1)),
        map=out[:, : out.shape[1] - extend].reshape(
            (1, out.shape[0], out.shape[1] - extend)
        ),
        spin=0,
        lmax=lmax,
        geometry="GL",
        nthreads=nthreads,
    )
    if extend:
        out[:, -1] = out[:, 0]
    return out


def _analyze_GLQ(map, lmax):
    alm = ducc0.sht.experimental.analysis_2d(
        map=map.reshape((1, map.shape[0], map.shape[1])),
        spin=0,
        lmax=lmax,
        geometry="GL",
        nthreads=nthreads,
    )
    return alm[0]


# ---- Real/Complex Coefficients transformations ----

def _ccilm2almr(cilm):
    lmax = cilm.shape[1] - 1
    alm = _np.empty((_nalm(lmax, lmax),), dtype=_np.complex128)
    fct = (-1) ** _np.arange(lmax + 1)
    alm[0: lmax + 1] = cilm[0, :, 0].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = _np.conj(cilm[1, m:, m])
        tmp *= fct[m]
        tmp += cilm[0, m:, m]
        tmp *= 1.0 / _np.sqrt(2.0)
        alm[ofs: ofs + lmax + 1 - m] = _np.conj(tmp)
        ofs += lmax + 1 - m
    return alm


def _ccilm2almi(cilm):
    lmax = cilm.shape[1] - 1
    alm = _np.empty((_nalm(lmax, lmax),), dtype=_np.complex128)
    fct = (-1) ** _np.arange(lmax + 1)
    alm[0: lmax + 1] = cilm[0, :, 0].imag
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = _np.conj(cilm[1, m:, m])
        tmp *= -fct[m]
        tmp += cilm[0, m:, m]
        tmp *= 1.0 / _np.sqrt(2.0)
        alm[ofs: ofs + lmax + 1 - m] = tmp.imag + 1j * tmp.real
        ofs += lmax + 1 - m
    return alm


def _addRealpart(cilm, alm):
    lmax = cilm.shape[1] - 1
    cilm[0, :, 0].real += alm[0: lmax + 1].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = alm[ofs: ofs + lmax + 1 - m] / _np.sqrt(2.0)
        cilm[0, m:, m].real += tmp.real
        cilm[0, m:, m].imag -= tmp.imag
        if m & 1:
            cilm[1, m:, m] -= tmp
        else:
            cilm[1, m:, m] += tmp
        ofs += lmax + 1 - m
    return cilm


def _addImagpart(cilm, alm):
    lmax = cilm.shape[1] - 1
    cilm[0, :, 0].imag += alm[0: lmax + 1].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = alm[ofs: ofs + lmax + 1 - m] / _np.sqrt(2.0)
        cilm[0, m:, m].real += tmp.imag
        cilm[0, m:, m].imag += tmp.real
        if m & 1:
            cilm[1, m:, m].real += tmp.imag
            cilm[1, m:, m].imag -= tmp.real
        else:
            cilm[1, m:, m].real -= tmp.imag
            cilm[1, m:, m].imag += tmp.real
        ofs += lmax + 1 - m
    return cilm


def _prep_lmax(lmax, lmax_calc, cilm):
    if lmax is None:
        lmax = cilm.shape[1] - 1
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    if lmax_calc > lmax:
        raise RuntimeError(
            "lmax_calc ({}) must be less than or equal to lmax ({})".format(
                lmax_calc, lmax
            )
        )
    # lmax_calc need must not be higher than cilm.shape[1] - 1.
    lmax_calc = min(cilm.shape[1] - 1, lmax_calc)
    return lmax, lmax_calc, cilm[:, : lmax_calc + 1, : lmax_calc + 1]


def SHRotateRealCoef(cilm, x, dj=None):
    """Determine the spherical harmonic coefficients of a real function rotated by
    three Euler angles.
    """ # noqa
    lmax = cilm.shape[1] - 1
    alm = _make_alm(cilm, lmax, 1, 1)
    alm = ducc0.sht.rotate_alm(alm, lmax, -x[0], -x[1], -x[2],
                               nthreads=nthreads)
    return _extract_alm(alm, lmax, 1, 1)


def SHRotateComplexCoef(cilm, x, dj=None):
    """Determine the spherical harmonic coefficients of a complex-valued function
    rotated by three Euler angles.
    """ # noqa
    lmax = cilm.shape[1] - 1
    alm = _ccilm2almr(cilm)
    alm = _apply_norm(alm, lmax, 1, 1, False)
    alm = ducc0.sht.rotate_alm(alm, lmax, -x[0], -x[1], -x[2],
                               nthreads=nthreads)
    alm = _apply_norm(alm, lmax, 1, 1, True)
    res = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.complex128)
    _addRealpart(res, alm)
    alm = _ccilm2almi(cilm)
    alm = _apply_norm(alm, lmax, 1, 1, False)
    alm = ducc0.sht.rotate_alm(alm, lmax, -x[0], -x[1], -x[2],
                               nthreads=nthreads)
    alm = _apply_norm(alm, lmax, 1, 1, True)
    _addImagpart(res, alm)
    return res


# ---- High-level DH functions ----

def MakeGridDH(
    cilm,
    lmax=None,
    norm=1,
    sampling=1,
    csphase=1,
    lmax_calc=None,
    extend=False,
):
    """Create a 2D map from a set of spherical harmonic coefficients using the Driscoll
    and Healy (1994) sampling theorem.
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _make_alm(cilm, lmax_calc, norm, csphase)
    out = _np.empty([2 * lmax + 2 + extend,
                     sampling * (2 * lmax + 2) + extend])
    return _synthesize_DH(alm, lmax_calc, extend, out)


def MakeGridDHC(
    cilm,
    lmax=None,
    norm=1,
    sampling=1,
    csphase=1,
    lmax_calc=None,
    extend=False,
):
    """Create a 2D complex map from a set of complex spherical harmonic coefficients
    that conforms with Driscoll and Healy's (1994) sampling theorem.
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _ccilm2almi(cilm)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, False)
    res = _np.empty(
        [2 * lmax + 2 + extend, sampling * (2 * lmax + 2) + extend],
        dtype=_np.complex128,
    )
    _synthesize_DH(alm, lmax_calc, extend, res.imag)
    alm = _ccilm2almr(cilm)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, False)
    _synthesize_DH(alm, lmax_calc, extend, res.real)
    return res


def SHExpandDH(griddh, norm=1, sampling=1, csphase=1, lmax_calc=None):
    """Expand an equally sampled or equally spaced grid into spherical harmonics using
    Driscoll and Healy's (1994) sampling theorem.
    """ # noqa
    griddh = _fixdtype(griddh)
    if griddh.shape[1] != sampling * griddh.shape[0]:
        raise RuntimeError("grid resolution mismatch")
    if lmax_calc is None:
        lmax_calc = griddh.shape[0] // 2 - 1
    if lmax_calc > (griddh.shape[0] // 2 - 1):
        raise RuntimeError("lmax_calc too high")
    alm = _analyze_DH(griddh, lmax_calc)
    return _extract_alm(alm, lmax_calc, norm, csphase)


def SHExpandDHC(griddh, norm=1, sampling=1, csphase=1, lmax_calc=None):
    """Expand an equally sampled or equally spaced complex grid into complex spherical
    harmonics using Driscoll and Healy's (1994) sampling theorem.
    """ # noqa
    if griddh.shape[1] != sampling * griddh.shape[0]:
        raise RuntimeError("grid resolution mismatch")
    if lmax_calc is None:
        lmax_calc = griddh.shape[0] // 2 - 1
    if lmax_calc > (griddh.shape[0] // 2 - 1):
        raise RuntimeError("lmax_calc too high")
    lmax = griddh.shape[0] // 2 - 1 if lmax_calc is None else lmax_calc
    res = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.complex128)
    alm = _analyze_DH(_fixdtype(griddh.real), lmax_calc)
    alm = _apply_norm(alm, lmax, norm, csphase, True)
    _addRealpart(res, alm)
    alm = _analyze_DH(_fixdtype(griddh.imag), lmax_calc)
    alm = _apply_norm(alm, lmax, norm, csphase, True)
    _addImagpart(res, alm)
    return res


def MakeGradientDH(
    cilm, lmax=None, sampling=1, lmax_calc=None, extend=False, radius=None
):
    """Compute the gradient of a scalar function and return grids of the two horizontal
    components that conform with Driscoll and Healy's (1994) sampling theorem.
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _make_alm(cilm, lmax_calc, 1, 1)
    res = _np.empty((2, 2 * lmax + 2 + extend,
                     sampling * (2 * lmax + 2) + extend))
    res = _synthesize_DH_deriv1(alm, lmax_calc, extend, res)
    if radius is not None:
        res *= 1.0 / radius
    return res


# ---- High-level CC functions ----

def MakeGridCC(
    cilm,
    lmax=None,
    norm=1,
    csphase=1,
    lmax_calc=None,
    extend=False,
):
    """Create a 2D map from a set of spherical harmonic coefficients using the
    Clenshaw-Curtis sampling theorem.
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _make_alm(cilm, lmax_calc, norm, csphase)
    ntheta = lmax + 2
    nphi = 2 * (lmax + 1)
    out = _np.empty([ntheta, nphi + extend])
    return _synthesize_CC(alm, lmax_calc, extend, out)


def MakeGridCCC(
    cilm,
    lmax=None,
    norm=1,
    csphase=1,
    lmax_calc=None,
    extend=False,
):
    """Create a 2D complex map from a set of complex spherical harmonic coefficients
    sampled on Clenshaw-Curtis nodes.
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _ccilm2almi(cilm)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, False)
    ntheta = lmax + 2
    nphi = 2 * (lmax + 1)
    res = _np.empty([ntheta, nphi + extend], dtype=_np.complex128)
    _synthesize_CC(alm, lmax_calc, extend, res.imag)
    alm = _ccilm2almr(cilm)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, False)
    _synthesize_CC(alm, lmax_calc, extend, res.real)
    return res


def SHExpandCC(gridcc, norm=1, csphase=1, lmax_calc=None):
    """Expand a Clenshaw-Curtis grid into spherical harmonics.
    """ # noqa
    gridcc = _fixdtype(gridcc)
    if lmax_calc is None:
        lmax_calc = gridcc.shape[0] - 2

    alm = _analyze_CC(gridcc, lmax_calc)
    return _extract_alm(alm, lmax_calc, norm, csphase)


def SHExpandCCC(gridcc, norm=1, csphase=1, lmax_calc=None):
    """Expand a complex Clenshaw-Curtis grid into complex spherical harmonics.
    """ # noqa
    if lmax_calc is None:
        lmax_calc = gridcc.shape[0] - 2

    lmax = gridcc.shape[0] - 2 if lmax_calc is None else lmax_calc
    res = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.complex128)
    alm = _analyze_CC(_fixdtype(gridcc.real), lmax_calc)
    alm = _apply_norm(alm, lmax, norm, csphase, True)
    _addRealpart(res, alm)
    alm = _analyze_CC(_fixdtype(gridcc.imag), lmax_calc)
    alm = _apply_norm(alm, lmax, norm, csphase, True)
    _addImagpart(res, alm)
    return res


def MakeGradientCC(
    cilm, lmax=None, lmax_calc=None, extend=False, radius=None
):
    """Compute the gradient of a scalar function and return grids of the two horizontal
    components that conform with Clenshaw-Curtis sampling.
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _make_alm(cilm, lmax_calc, 1, 1)
    ntheta = lmax + 2
    nphi = 2 * (lmax + 1)
    res = _np.empty((2, ntheta, nphi + extend))
    res = _synthesize_CC_deriv1(alm, lmax_calc, extend, res)
    if radius is not None:
        res *= 1.0 / radius
    return res


def MakeGradientCCC(
    cilm,
    lmax=None,
    norm=1,
    csphase=1,
    lmax_calc=None,
    extend=False,
    radius=None
):
    """Compute the gradient of a complex scalar function and return grids of the
    two horizontal complex components that conform with Clenshaw-Curtis sampling.
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    ntheta = lmax + 2
    nphi = 2 * (lmax + 1)
    res = _np.empty((2, ntheta, nphi + extend), dtype=_np.complex128)

    # Real part of the function
    alm = _ccilm2almr(cilm)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, False)
    tmp = _np.empty((2, ntheta, nphi + extend))
    _synthesize_CC_deriv1(alm, lmax_calc, extend, tmp)
    res.real = tmp

    # Imaginary part of the function
    alm = _ccilm2almi(cilm)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, False)
    _synthesize_CC_deriv1(alm, lmax_calc, extend, tmp)
    res.imag = tmp

    if radius is not None:
        res *= 1.0 / radius
    return res


# ---- High-level GLQ functions ----

def MakeGridGLQ(
    cilm, zero=None, lmax=None, norm=1, csphase=1, lmax_calc=None, extend=False
):
    """Create a 2D map from a set of spherical harmonic coefficients sampled on the
    Gauss-Legendre quadrature nodes.
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _make_alm(cilm, lmax_calc, norm, csphase)
    out = _np.empty([lmax + 1, (2 * lmax + 1) + extend])
    return _synthesize_GLQ(alm, lmax_calc, extend, out)


def MakeGridGLQC(
    cilm, zero=None, lmax=None, norm=1, csphase=1, lmax_calc=None, extend=False
):
    """Create a 2D complex map from a set of complex spherical harmonic coefficients
    sampled on the Gauss-Legendre quadrature nodes.
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _ccilm2almi(cilm)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, False)
    res = _np.empty([lmax + 1, 2 * lmax + 1 + extend], dtype=_np.complex128)
    _synthesize_GLQ(alm, lmax_calc, extend, res.imag)
    alm = _ccilm2almr(cilm)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, False)
    _synthesize_GLQ(alm, lmax_calc, extend, res.real)
    return res


def SHExpandGLQ(gridglq, w=None, zero=None, norm=1, csphase=1, lmax_calc=None):
    """
    Expand a 2D grid sampled on the Gauss-Legendre quadrature nodes into spherical
    harmonics.
    """ # noqa
    gridglq = _fixdtype(gridglq)
    if lmax_calc is None:
        lmax_calc = gridglq.shape[0] - 1
    if lmax_calc > (gridglq.shape[0] - 1):
        raise RuntimeError("lmax_calc too high")
    alm = _analyze_GLQ(gridglq, lmax_calc)
    return _extract_alm(alm, lmax_calc, norm, csphase)


def SHExpandGLQC(gridglq, w=None, zero=None, norm=1, csphase=1,
                 lmax_calc=None):
    """Expand a 2D grid sampled on the Gauss-Legendre quadrature nodes into spherical
    harmonics.
    """ # noqa
    if lmax_calc is None:
        lmax_calc = gridglq.shape[0] - 1
    if lmax_calc > (gridglq.shape[0] - 1):
        raise RuntimeError("lmax_calc too high")
    res = _np.zeros((2, lmax_calc + 1, lmax_calc + 1), dtype=_np.complex128)
    alm = _analyze_GLQ(_fixdtype(gridglq.real), lmax_calc)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, True)
    _addRealpart(res, alm)
    alm = _analyze_GLQ(_fixdtype(gridglq.imag), lmax_calc)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, True)
    _addImagpart(res, alm)
    return res
