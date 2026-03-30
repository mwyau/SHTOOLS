"""
This script tests the vector spherical harmonics.
"""
import numpy as np
import pyshtools as pysh

pysh.utils.figstyle()

# ==== R_EARTH, LMAX, OMEGA ====
R_EARTH = pysh.constants.Earth.wgs84.r3.value
LMAX = 71
OMEGA = 2 * np.pi / (24 * 60 * 60)


def main():
    test_solid_body_rotation()
    test_pure_divergence()
    test_round_trip()
    test_harmonic_round_trip()
    test_gradient_streamfunction()
    test_backend_enforcement()


def test_solid_body_rotation():
    """
    Test solid body rotation.
    solid body rotation, u = R * omega * sin(colat), v = 0
    """
    lmax = LMAX
    grid = pysh.SHGrid.from_zeros(lmax, grid='CC')
    lats = grid.lats()
    lons = grid.lons()
    colat_1d = np.radians(90. - lats)

    # Broadcast to 2D
    u_1d = R_EARTH * OMEGA * np.sin(colat_1d)
    u = np.broadcast_to(u_1d[:, np.newaxis], (len(lats), len(lons)))
    v = np.zeros_like(u)
    vec = pysh.SHVectorGrid.from_uv(u, v, a=R_EARTH)

    coeffs = vec.expand()
    div = coeffs.divergence().expand(grid='CC')
    vort = coeffs.vorticity().expand(grid='CC')

    # divergence should be zero
    np.testing.assert_allclose(div.data, 0, atol=1e-12)

    # vorticity should be -2 * omega * cos(colat)
    analytical_vort_1d = -2 * OMEGA * np.cos(colat_1d)
    analytical_vort = np.broadcast_to(analytical_vort_1d[:, np.newaxis],
                                      vort.data.shape)
    np.testing.assert_allclose(vort.data, analytical_vort, atol=1e-12)


def test_pure_divergence():
    """
    Test pure divergence.
    """
    lmax = LMAX
    grid = pysh.SHGrid.from_zeros(lmax, grid='CC')
    lats = grid.lats()
    lons = grid.lons()
    colat_1d = np.radians(90. - lats)

    v_1d = np.sin(colat_1d)
    v = np.broadcast_to(v_1d[:, np.newaxis], (len(lats), len(lons)))
    u = np.zeros_like(v)
    vec = pysh.SHVectorGrid.from_uv(u, v, a=R_EARTH)

    coeffs = vec.expand()
    div = coeffs.divergence().expand(grid='CC')
    vort = coeffs.vorticity().expand(grid='CC')

    # vorticity should be zero
    np.testing.assert_allclose(vort.data, 0, atol=1e-12)

    # divergence should be 2 / R * cos(colat)
    analytical_div_1d = 2.0 / R_EARTH * np.cos(colat_1d)
    analytical_div = np.broadcast_to(analytical_div_1d[:, np.newaxis],
                                     div.data.shape)
    np.testing.assert_allclose(div.data, analytical_div, atol=1e-12)


def test_round_trip():
    """
    Test round-trip reversibility with smooth field.
    """
    lmax = LMAX
    # Create a smooth field by starting from coefficients
    power = np.ones(lmax + 1)
    coeffs_e = pysh.SHCoeffs.from_random(power)
    coeffs_b = pysh.SHCoeffs.from_random(power)
    # Taper coefficients to make it very smooth
    l = np.arange(lmax + 1)
    taper = np.exp(-l**2 / (2 * 10**2))
    for l_val in range(lmax + 1):
        coeffs_e.coeffs[:, l_val, :] *= taper[l_val]
        coeffs_b.coeffs[:, l_val, :] *= taper[l_val]

    vec_coeffs = pysh.SHVectorCoeffs(coeffs_e, coeffs_b, a=R_EARTH)
    vec_grid = vec_coeffs.expand(grid='CC')
    vec_coeffs_back = vec_grid.expand()

    np.testing.assert_allclose(vec_coeffs.e_mode.coeffs[:, 1:, :],
                               vec_coeffs_back.e_mode.coeffs[:, 1:, :],
                               atol=1e-12)
    np.testing.assert_allclose(vec_coeffs.b_mode.coeffs[:, 1:, :],
                               vec_coeffs_back.b_mode.coeffs[:, 1:, :],
                               atol=1e-12)


def test_harmonic_round_trip():
    """
    Test round-trip for a single harmonic.
    """
    lmax = LMAX
    coeffs_e = pysh.SHCoeffs.from_zeros(lmax)
    coeffs_e.set_coeffs(1.0, 2, 1)
    coeffs_b = pysh.SHCoeffs.from_zeros(lmax)

    vec_coeffs = pysh.SHVectorCoeffs(coeffs_e, coeffs_b, a=R_EARTH)
    vec_grid = vec_coeffs.expand(grid='CC')
    vec_coeffs_back = vec_grid.expand()

    np.testing.assert_allclose(vec_coeffs.e_mode.coeffs,
                               vec_coeffs_back.e_mode.coeffs, atol=1e-12)


def test_gradient_streamfunction():
    """
    Test gradient_vector and streamfunction_vector methods.
    """
    lmax = LMAX
    coeffs = pysh.SHCoeffs.from_zeros(lmax)
    coeffs.set_coeffs(1.0, 2, 0)  # Y_20 = 1/4 sqrt(5/pi) (3 cos^2 theta - 1)

    # Gradient
    v_grad = coeffs.gradient_vector(radius=R_EARTH)
    div = v_grad.divergence()
    # div should be Laplacian(Y_20) / R^2 = -l(l+1)/R^2 Y_20 = -6/R^2 Y_20
    expected_div = coeffs.copy()
    expected_div.coeffs *= -6.0 / R_EARTH**2
    np.testing.assert_allclose(div.coeffs, expected_div.coeffs, atol=1e-12)

    # Streamfunction
    v_stream = coeffs.streamfunction_vector(radius=R_EARTH)
    vort = v_stream.vorticity()
    # vort should be Laplacian(Y_20) / R^2 = -6/R^2 Y_20
    np.testing.assert_allclose(vort.coeffs, expected_div.coeffs, atol=1e-12)


def test_backend_enforcement():
    """
    Test backend enforcement.
    """
    pysh.backends.select_preferred_backend('shtools')
    lmax = 1
    grid = pysh.SHGrid.from_zeros(lmax, grid='CC')
    u = np.zeros((grid.nlat, grid.nlon))
    v = np.zeros((grid.nlat, grid.nlon))
    vec = pysh.SHVectorGrid.from_uv(u, v, a=R_EARTH)

    try:
        vec.expand()
    except NotImplementedError:
        pass
    else:
        raise RuntimeError("Should have raised NotImplementedError for "
                           "shtools backend")

    pysh.backends.select_preferred_backend('ducc')


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
