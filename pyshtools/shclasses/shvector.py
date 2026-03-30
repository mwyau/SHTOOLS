"""
Vector spherical harmonic classes.
"""
import numpy as np
import pyshtools as pysh
import copy as _copy
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import xarray as _xr

from ..backends import ducc0_wrapper as ducc
from .shgrid import SHGrid as _SHGrid
from .shcoeffs import SHCoeffs as _SHCoeffs
from ..constants.Earth import wgs84 as _wgs84

try:
    import cartopy.crs as _ccrs
    _cartopy_module = True
except ImportError:
    _cartopy_module = False


class SHVectorGrid(object):
    """
    Class for vector fields on the sphere.

    The vector field is represented by its meridional (vtheta) and
    zonal (vphi) components. In terms of zonal (u) and meridional (v)
    components, vtheta = v and vphi = u.

    Attributes:

    vtheta         : SHGrid class instance of the meridional component.
    vphi           : SHGrid class instance of the zonal component.
    a              : The radius of the sphere.
    lmax           : The maximum spherical harmonic degree resolvable by the
                     grids.
    nlat, nlon     : The number of latitude and longitude bands in the grids.
    grid           : The grid type ('DH', 'GLQ', or 'CC').

    Methods:

    u()            : Return the zonal component of the vector field.
    v()            : Return the meridional component of the vector field.
    expand()       : Expand the vector field into spherical harmonic
                     coefficients.
    plot()         : Plot the two components of the vector field.
    plot_vtheta()  : Plot the meridional component of the vector field.
    plot_vphi()    : Plot the zonal component of the vector field.
    quiver()       : Plot the vector field as quivers.
    to_xarray()    : Return an xarray DataSet of the vector field.
    copy()         : Return a copy of the class instance.
    info()         : Print a summary of the data stored in the SHVectorGrid
                     instance.
    """
    def __init__(self, vtheta, vphi, a=None, grid='CC'):
        """
        Initialize the vector class.
        """
        if grid.upper() not in ['DH', 'GLQ', 'CC']:
            raise ValueError(
                "grid must be 'DH', 'GLQ' or 'CC'. " +
                "Input value is {:s}".format(repr(grid)))

        self.grid = grid.upper()

        if isinstance(vtheta, _SHGrid):
            self.vtheta = vtheta
        else:
            self.vtheta = _SHGrid.from_array(vtheta, grid=self.grid)

        if isinstance(vphi, _SHGrid):
            self.vphi = vphi
        else:
            self.vphi = _SHGrid.from_array(vphi, grid=self.grid)

        if a is None:
            self.a = _wgs84.r3.value
        else:
            self.a = a

        self.lmax = self.vtheta.lmax
        self.nlat = self.vtheta.nlat
        self.nlon = self.vtheta.nlon

    @classmethod
    def from_uv(cls, u, v, a=None, grid='CC'):
        """
        Initialize the vector class from the zonal and meridional components.

        Note: ducc0 spin-1 expects (v_theta, v_phi). For meteorological
        conventions where u is zonal and v is northward, v_theta=v and
        v_phi=u matches the convention used in NCL's Spherepack wrappers.
        """
        vtheta = v
        vphi = u
        return cls(vtheta, vphi, a=a, grid=grid)

    def u(self):
        """
        Return the zonal component of the vector field.
        """
        return self.vphi

    def v(self):
        """
        Return the meridional (northward) component of the vector field.
        """
        return self.vtheta

    def copy(self):
        """
        Return a deep copy of the class instance.
        """
        return _copy.deepcopy(self)

    def info(self):
        """
        Print a summary of the data stored in the SHVectorGrid instance.
        """
        print(repr(self))

    def __repr__(self):
        return (f'  grid = {self.grid!r}\n'
                f'  nlat = {self.nlat}\n'
                f'  nlon = {self.nlon}\n'
                f'  lmax = {self.lmax}\n'
                f'  radius = {self.a}\n'
                )

    def plot_vtheta(self, **kwargs):
        """
        Plot the meridional component of the vector field.
        """
        if 'cb_label' not in kwargs:
            kwargs['cb_label'] = r'$v_\theta$'
        return self.vtheta.plot(**kwargs)

    def plot_vphi(self, **kwargs):
        """
        Plot the zonal component of the vector field.
        """
        if 'cb_label' not in kwargs:
            kwargs['cb_label'] = r'$v_\phi$'
        return self.vphi.plot(**kwargs)

    def plot(self, colorbar='bottom', **kwargs):
        """
        Plot the two components of the vector field.
        """
        if colorbar is not None:
            if colorbar in set(['bottom', 'top']):
                scale = 0.4
            else:
                scale = 0.25
        else:
            scale = 0.3
        figsize = (_mpl.rcParams['figure.figsize'][0],
                   _mpl.rcParams['figure.figsize'][0] * scale)

        fig, ax = _plt.subplots(1, 2, figsize=figsize)

        kwargs_theta = kwargs.copy()
        if 'cb_label' not in kwargs_theta:
            kwargs_theta['cb_label'] = r'$v_\theta$'
        self.vtheta.plot(ax=ax[0], colorbar=colorbar, show=False,
                         **kwargs_theta)

        kwargs_phi = kwargs.copy()
        if 'cb_label' not in kwargs_phi:
            kwargs_phi['cb_label'] = r'$v_\phi$'
        self.vphi.plot(ax=ax[1], colorbar=colorbar, show=False,
                       **kwargs_phi)

        fig.tight_layout(pad=0.5)
        if kwargs.get('show', True):
            _plt.show()
        return fig, ax

    def quiver(self, ax=None, projection=None, step=1, **kwargs):
        """
        Plot the vector field as quivers.

        Usage
        -----
        ax = x.quiver([ax, projection, step, **kwargs])

        Parameters
        ----------
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        projection : Cartopy projection class, optional, default = None
            The Cartopy projection class used to project the data.
        step : int, optional, default = 1
            The step size for subsampling the grid.
        **kwargs : optional
            Optional arguments passed to matplotlib.pyplot.quiver().
        """
        if ax is None:
            if projection is not None:
                if not _cartopy_module:
                    raise ImportError('quiver() requires Cartopy for '
                                      'map projections.')
                fig, ax = _plt.subplots(subplot_kw={'projection': projection})
            else:
                fig, ax = _plt.subplots()

        lats = self.vtheta.lats()
        lons = self.vphi.lons()
        lon2d, lat2d = np.meshgrid(lons, lats)

        sl = slice(None, None, step)
        u = self.u().data[sl, sl]
        v = self.v().data[sl, sl]

        if projection is not None:
            ax.quiver(lon2d[sl, sl], lat2d[sl, sl], u, v,
                      transform=_ccrs.PlateCarree(), **kwargs)
        else:
            ax.quiver(lon2d[sl, sl], lat2d[sl, sl], u, v, **kwargs)

        return ax

    def to_xarray(self, title='', description='',
                  comment='pyshtools vector grid'):
        """
        Return the vector field gridded data as an xarray DataSet.
        """
        attrs = {'title': title,
                 'description': description,
                 'comment': comment,
                 'nlat': self.nlat,
                 'nlon': self.nlon,
                 'lmax': self.lmax,
                 'grid': self.grid,
                 'radius': self.a
                 }

        _vtheta = self.vtheta.to_xarray(title='v_theta',
                                        long_name='meridional component')
        _vphi = self.vphi.to_xarray(title='v_phi',
                                    long_name='zonal component')

        return _xr.Dataset({'vtheta': _vtheta, 'vphi': _vphi}, attrs=attrs)

    def expand(self, lmax_calc=None):
        """
        Expand the vector field into its spherical harmonic coefficients.
        """
        if pysh.backends.preferred_backend() != 'ducc':
            raise NotImplementedError(
                "Vector transforms are only implemented for the 'ducc' " +
                "backend.")

        lmax = self.lmax if lmax_calc is None else lmax_calc

        map_spin1 = np.stack((self.vtheta.data, self.vphi.data), axis=0)

        # Dynamic detection of longitudinal extension:
        if self.grid in ('DH', 'CC'):
            if map_spin1.shape[2] % 2 != 0:
                map_spin1 = map_spin1[:, :, :-1]
        else:  # GLQ
            if map_spin1.shape[2] % 2 == 0:
                map_spin1 = map_spin1[:, :, :-1]

        if self.grid == 'DH':
            geometry = 'DH'
        elif self.grid == 'CC':
            geometry = 'CC'
        else:
            geometry = 'GL'

        alm_E, alm_B = ducc.ducc0.sht.experimental.analysis_2d(
            map=map_spin1, spin=1, lmax=lmax, geometry=geometry,
            nthreads=ducc.nthreads)

        cilm_e = ducc._extract_alm(alm_E, lmax, norm=1, csphase=1)
        cilm_b = ducc._extract_alm(alm_B, lmax, norm=1, csphase=1)

        return SHVectorCoeffs(cilm_e, cilm_b, a=self.a)


class SHVectorCoeffs(object):
    """
    Class for vector spherical harmonic coefficients.

    The vector field is represented by its E-mode (divergence-like) and
    B-mode (vorticity-like) coefficients.

    Attributes:

    e_mode         : SHCoeffs class instance of the E-mode coefficients.
    b_mode         : SHCoeffs class instance of the B-mode coefficients.
    a              : The radius of the sphere.
    lmax           : The maximum spherical harmonic degree of the coefficients.

    Methods:

    from_array()   : Initialize the vector class from two arrays of spherical
                     harmonic coefficients.
    from_vorticity_divergence() : Initialize from vorticity and divergence.
    from_potential_streamfunction() : Initialize from velocity potential and
                                      streamfunction.
    divergence()   : Return the divergence of the vector field.
    vorticity()    : Return the vorticity of the vector field.
    velocity_potential() : Return the velocity potential.
    streamfunction()     : Return the streamfunction.
    expand()       : Expand the coefficients into a vector field.
    copy()         : Return a copy of the class instance.
    info()         : Print a summary of the data stored in the SHVectorCoeffs
                     instance.
    """
    def __init__(self, coeffs_e, coeffs_b, a=None):
        """
        Initialize the vector class from spherical harmonic coefficients.
        """
        if hasattr(coeffs_e, 'coeffs'):
            self.e_mode = coeffs_e
        else:
            self.e_mode = _SHCoeffs.from_array(coeffs_e)

        if hasattr(coeffs_b, 'coeffs'):
            self.b_mode = coeffs_b
        else:
            self.b_mode = _SHCoeffs.from_array(coeffs_b)

        if a is None:
            self.a = _wgs84.r3.value
        else:
            self.a = a

        self.lmax = self.e_mode.lmax

    @classmethod
    def from_array(cls, coeffs_e, coeffs_b, a=None, **kwargs):
        """
        Initialize the vector class from two arrays of spherical harmonic
        coefficients.
        """
        e_mode = _SHCoeffs.from_array(coeffs_e, **kwargs)
        b_mode = _SHCoeffs.from_array(coeffs_b, **kwargs)
        return cls(e_mode, b_mode, a=a)

    @classmethod
    def from_vorticity_divergence(cls, vorticity, divergence, a=None):
        """
        Initialize the vector class from vorticity and divergence.
        """
        if isinstance(vorticity, _SHCoeffs):
            vort = vorticity
        else:
            vort = _SHCoeffs.from_array(vorticity)

        if isinstance(divergence, _SHCoeffs):
            div = divergence
        else:
            div = _SHCoeffs.from_array(divergence)

        if a is None:
            a = _wgs84.r3.value

        lmax = vort.lmax
        l = np.arange(lmax + 1)
        # Avoid division by zero for l=0
        inv_scale = np.zeros(lmax + 1)
        inv_scale[1:] = a / np.sqrt(l[1:] * (l[1:] + 1))

        # Div = -sqrt(l(l+1))/R * E => E = -Div * R / sqrt(l(l+1))
        # Vort = -sqrt(l(l+1))/R * B => B = -Vort * R / sqrt(l(l+1))
        coeffs_e = div.copy()
        coeffs_e.coeffs *= -1
        coeffs_b = vort.copy()
        coeffs_b.coeffs *= -1

        for l_val in range(1, lmax + 1):
            coeffs_e.coeffs[:, l_val, :l_val+1] *= inv_scale[l_val]
            coeffs_b.coeffs[:, l_val, :l_val+1] *= inv_scale[l_val]

        return cls(coeffs_e, coeffs_b, a=a)

    @classmethod
    def from_potential_streamfunction(cls, potential, streamfunction, a=None):
        """
        Initialize the vector class from velocity potential and streamfunction.
        """
        if isinstance(potential, _SHCoeffs):
            pot = potential
        else:
            pot = _SHCoeffs.from_array(potential)

        if isinstance(streamfunction, _SHCoeffs):
            stream = streamfunction
        else:
            stream = _SHCoeffs.from_array(streamfunction)

        if a is None:
            a = _wgs84.r3.value

        lmax = pot.lmax
        l = np.arange(lmax + 1)
        scale = np.sqrt(l * (l + 1)) / a

        # Chi = -R/sqrt(l(l+1)) * E => E = -sqrt(l(l+1))/R * Chi
        # Psi = -R/sqrt(l(l+1)) * B => B = -sqrt(l(l+1))/R * Psi
        coeffs_e = pot.copy()
        coeffs_e.coeffs *= -1
        coeffs_b = stream.copy()
        coeffs_b.coeffs *= -1

        for l_val in range(1, lmax + 1):
            coeffs_e.coeffs[:, l_val, :l_val+1] *= scale[l_val]
            coeffs_b.coeffs[:, l_val, :l_val+1] *= scale[l_val]

        return cls(coeffs_e, coeffs_b, a=a)

    def copy(self):
        """
        Return a deep copy of the class instance.
        """
        return _copy.deepcopy(self)

    def info(self):
        """
        Print a summary of the data stored in the SHVectorCoeffs instance.
        """
        print(repr(self))

    def __repr__(self):
        return (f'  lmax = {self.lmax}\n'
                f'  radius = {self.a}\n'
                f'  normalization = {self.e_mode.normalization!r}\n'
                f'  csphase = {self.e_mode.csphase}\n'
                )

    def _vorticity_divergence(self, kind):
        """
        Return the vorticity or divergence coefficients of the vector field.
        """
        if kind == 'vorticity':
            # Vorticity = sqrt(l(l+1))/R * B = Laplacian(psi)
            coeffs = self.b_mode.copy()
        else:
            # Divergence = -sqrt(l(l+1))/R * E = Laplacian(chi)
            coeffs = self.e_mode.copy()
            coeffs.coeffs *= -1

        l = np.arange(self.lmax + 1)
        scale = np.sqrt(l * (l + 1)) / self.a
        scale[0] = 0.
        scaled_coeffs = coeffs.copy()
        for l_val in range(1, self.lmax + 1):
            scaled_coeffs.coeffs[:, l_val, :l_val+1] *= scale[l_val]

        return scaled_coeffs

        return scaled_coeffs

    def divergence(self):
        """
        Return the divergence of the vector field.
        """
        return self._vorticity_divergence('divergence')

    def vorticity(self):
        """
        Return the vorticity of the vector field.
        """
        return self._vorticity_divergence('vorticity')

    def _potential(self, kind):
        """
        Return the velocity potential or streamfunction coefficients.
        """
        if kind == 'potential':
            # Chi = -R/sqrt(l(l+1)) * E
            coeffs = self.e_mode.copy()
        else:
            # Psi = -R/sqrt(l(l+1)) * B
            coeffs = self.b_mode.copy()
            coeffs.coeffs *= -1

        l = np.arange(self.lmax + 1)
        # Avoid division by zero for l=0
        inv_scale = np.zeros(self.lmax + 1)
        inv_scale[1:] = self.a / np.sqrt(l[1:] * (l[1:] + 1))

        scaled_coeffs = coeffs.copy()
        for l_val in range(1, self.lmax + 1):
            scaled_coeffs.coeffs[:, l_val, :l_val+1] *= inv_scale[l_val]

        return scaled_coeffs

    def velocity_potential(self):
        """
        Return the velocity potential coefficients of the vector field.
        """
        return self._potential('potential')

    def streamfunction(self):
        """
        Return the streamfunction coefficients of the vector field.
        """
        return self._potential('streamfunction')

    def expand(self, grid='CC', lmax_calc=None, extend=True):
        """
        Expand the spherical harmonic coefficients into a vector field.
        """
        if pysh.backends.preferred_backend() != 'ducc':
            raise NotImplementedError(
                "Vector transforms are only implemented for the 'ducc' " +
                "backend.")

        if grid.upper() not in ['DH', 'GLQ', 'CC']:
            raise ValueError(
                "grid must be 'DH', 'GLQ' or 'CC'. " +
                "Input value is {:s}".format(repr(grid)))

        lmax = self.lmax if lmax_calc is None else lmax_calc

        # ducc0 expects normalization=1 (4pi) and csphase=1
        coeffs_e_4pi = self.e_mode.to_array(normalization='4pi', csphase=1)
        coeffs_b_4pi = self.b_mode.to_array(normalization='4pi', csphase=1)

        alm_E = ducc._make_alm(coeffs_e_4pi, lmax, norm=1, csphase=1)
        alm_B = ducc._make_alm(coeffs_b_4pi, lmax, norm=1, csphase=1)
        alm_spin1 = np.stack((alm_E, alm_B), axis=0)

        if grid.upper() == 'DH':
            ntheta = 2 * lmax + 2
            nphi = 2 * ntheta
            geometry = 'DH'
            if extend:
                ntheta += 1
                nphi += 1
        elif grid.upper() == 'CC':
            ntheta = lmax + 2
            nphi = 2 * (lmax + 1)
            geometry = 'CC'
            if extend:
                nphi += 1
        else:  # GLQ
            ntheta = lmax + 1
            nphi = 2 * lmax + 1
            geometry = 'GL'
            if extend:
                nphi += 1

        vtheta_ducc, vphi_ducc = ducc.ducc0.sht.experimental.synthesis_2d(
            alm=alm_spin1.reshape((2, -1)), spin=1, lmax=lmax,
            ntheta=ntheta, nphi=nphi - (1 if extend else 0),
            geometry=geometry, nthreads=ducc.nthreads)

        if extend:
            vtheta = np.zeros((ntheta, nphi))
            vtheta[:, :-1] = vtheta_ducc
            vtheta[:, -1] = vtheta_ducc[:, 0]

            vphi = np.zeros((ntheta, nphi))
            vphi[:, :-1] = vphi_ducc
            vphi[:, -1] = vphi_ducc[:, 0]
        else:
            vtheta = vtheta_ducc
            vphi = vphi_ducc

        return SHVectorGrid(vtheta, vphi, a=self.a, grid=grid)
