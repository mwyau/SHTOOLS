# MakeGridCC()

Create a real Clenshaw-Curtis grid from spherical harmonic coefficients.

# Usage

gridcc = MakeGridCC(cilm, [lmax, norm, csphase, lmax_calc, extend])

# Returns

gridcc : float, dimension (lmax+2, 2\*lmax+2+extend)
:   A Clenshaw-Curtis grid containing both poles. If `extend` is True, the
    final column is the redundant 360 E longitude.

# Parameters

cilm : float, dimension (2, lmaxin+1, lmaxin+1)
:   Real cosine and sine spherical harmonic coefficients.

lmax : optional, integer, default = lmaxin
:   The maximum degree defining the output grid.

norm : optional, integer, default = 1
:   1 = geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized
    harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

csphase : optional, integer, default = 1
:   1 = exclude the Condon-Shortley phase factor; -1 = include it.

lmax_calc : optional, integer, default = lmax
:   The maximum degree synthesized. It must be between 0 and lmax.

extend : optional, bool, default = False
:   If True, append the redundant longitude column at 360 E.

# Description

`MakeGridCC` is the real inverse of `SHExpandCC`. The non-extended grid has
shape (lmax+2, 2\*lmax+2), includes both poles, and has no `sampling`
parameter. It evaluates the associated Legendre basis at equally spaced
colatitudes and uses an FFT in longitude. Unnormalized associated Legendre
functions are accurate only to about degree 15.
