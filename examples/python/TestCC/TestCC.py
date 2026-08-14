"""Regression tests for Clenshaw-Curtis scalar grids."""

import numpy as np
import pyshtools as pysh

from pyshtools.backends import ducc0_wrapper


def main():
    test_grid()
    test_real_transform()
    test_complex_transform()
    test_full_bandwidth()
    test_conventions()
    test_lmax_calc()
    test_cross_backend()


def test_grid():
    """Check CC grid dimensions, coordinates, poles, and extension."""
    lmax = 5
    for extend in (False, True):
        grid = pysh.SHGrid.from_zeros(lmax, grid='CC', extend=extend)
        expected_nlon = 2 * lmax + 2 + extend
        if grid.grid != 'CC' or grid.data.shape != (lmax + 2, expected_nlon):
            raise Exception('CC grid dimensions or grid type are incorrect.')
        if hasattr(grid, 'sampling'):
            raise Exception('CC grids must not define DH sampling.')
        np.testing.assert_allclose(grid.lats(),
                                   np.linspace(90.0, -90.0, lmax + 2))
        np.testing.assert_allclose(grid.lons(), np.linspace(
            0.0, 360.0 if extend else 330.0, expected_nlon))
        if extend:
            np.testing.assert_allclose(grid.lons()[-1], 360.0)
        else:
            np.testing.assert_allclose(grid.lons()[-1], 330.0)

        restored = pysh.SHGrid.from_xarray(grid.to_xarray(), grid='CC')
        if restored.grid != 'CC' or restored.extend != extend:
            raise Exception(
                'CC xarray construction did not preserve the grid.')
        np.testing.assert_allclose(restored.data, grid.data)


def test_real_transform():
    """Round-trip real coefficients using the native CC implementation."""
    coeffs = pysh.SHCoeffs.from_array(_real_coefficients(8))
    for extend in (False, True):
        grid = coeffs.expand(grid='CC', backend='shtools', extend=extend)
        if extend:
            np.testing.assert_allclose(grid.data[:, 0], grid.data[:, -1])
        _assert_coefficients(grid.expand(backend='shtools'), coeffs)


def test_complex_transform():
    """Round-trip complex coefficients using the native CC implementation."""
    coeffs = pysh.SHCoeffs.from_array(_complex_coefficients(8))
    grid = coeffs.expand(grid='CC', backend='shtools', extend=False)
    _assert_coefficients(grid.expand(backend='shtools'), coeffs)


def test_full_bandwidth():
    """Retain sectoral and near-sectoral modes that detect a DH fallback."""
    lmax = 8
    for complex_ in (False, True):
        coefficients = np.zeros(
            (2, lmax + 1, lmax + 1),
            dtype=np.complex128 if complex_ else np.float64,
        )
        coefficients[0, lmax, lmax] = 1.0 + 0.25j if complex_ else 1.0
        coefficients[1, lmax, lmax - 1] = -0.5 + 0.75j if complex_ else -0.5
        coeffs = pysh.SHCoeffs.from_array(coefficients)
        grid = coeffs.expand(grid='CC', backend='shtools', extend=False)
        _assert_coefficients(grid.expand(backend='shtools'), coeffs,
                             atol=3e-10)


def test_conventions():
    """Exercise CC round trips for normalization and phase conventions."""
    coefficients = _real_coefficients(6)
    for normalization in ('4pi', 'schmidt', 'unnorm', 'ortho'):
        atol = 2e-8 if normalization == 'unnorm' else 2e-10
        for csphase in (-1, 1):
            coeffs = pysh.SHCoeffs.from_array(
                coefficients, normalization=normalization, csphase=csphase)
            grid = coeffs.expand(grid='CC', backend='shtools', extend=False)
            _assert_coefficients(grid.expand(backend='shtools'), coeffs,
                                 atol=atol)


def test_lmax_calc():
    """Check that a reduced calculation degree returns the requested band."""
    coefficients = _real_coefficients(8)
    coeffs = pysh.SHCoeffs.from_array(coefficients)
    grid = coeffs.expand(grid='CC', lmax_calc=3, backend='shtools',
                         extend=False)
    recovered = grid.expand(lmax_calc=3, backend='shtools')
    np.testing.assert_allclose(
        recovered.to_array(normalization='4pi', csphase=1),
        coefficients[:, :4, :4], atol=2e-11)


def test_cross_backend():
    """Compare the native and DUCC CC paths when DUCC is available."""
    if not ducc0_wrapper.available():
        print('DUCC is unavailable; skipping CC cross-backend comparison.')
        return

    for coefficients in (_real_coefficients(8), _complex_coefficients(8)):
        coeffs = pysh.SHCoeffs.from_array(coefficients)
        native = coeffs.expand(grid='CC', backend='shtools', extend=False)
        ducc = coeffs.expand(grid='CC', backend='ducc', extend=False)
        np.testing.assert_allclose(ducc.data, native.data, atol=3e-10)
        _assert_coefficients(native.expand(backend='ducc'), coeffs,
                             atol=3e-10)


def _real_coefficients(lmax, seed=4):
    rng = np.random.default_rng(seed)
    coefficients = rng.standard_normal((2, lmax + 1, lmax + 1))
    for degree in range(lmax + 1):
        coefficients[:, degree, degree + 1:] = 0.0
    coefficients[1, :, 0] = 0.0
    return coefficients


def _complex_coefficients(lmax, seed=5):
    rng = np.random.default_rng(seed)
    coefficients = _real_coefficients(lmax, seed=seed).astype(np.complex128)
    coefficients += 1j * rng.standard_normal(coefficients.shape)
    for degree in range(lmax + 1):
        coefficients[:, degree, degree + 1:] = 0.0
    coefficients[1, :, 0] = 0.0
    return coefficients


def _assert_coefficients(recovered, expected, atol=2e-10):
    np.testing.assert_allclose(
        recovered.to_array(normalization=expected.normalization,
                           csphase=expected.csphase),
        expected.to_array(normalization=expected.normalization,
                          csphase=expected.csphase), atol=atol)


if __name__ == '__main__':
    main()
