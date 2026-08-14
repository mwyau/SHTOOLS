# SHExpandCCC()

Expand a complex Clenshaw-Curtis grid into spherical harmonics.

# Usage

cilm = SHExpandCCC(gridcc, [norm, csphase, lmax_calc])

# Returns

cilm : complex, dimension (2, lmax_calc+1, lmax_calc+1)
:   Complex coefficients following the `SHExpandDHC` convention: positive
    orders in the first index and negative orders in the second before the
    conventional `(-1)^m` factor.

# Parameters

gridcc : complex, dimension (L+2, 2\*L+2)
:   A non-extended Clenshaw-Curtis grid containing both poles.

norm : optional, integer, default = 1
:   1 = geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized
    harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

csphase : optional, integer, default = 1
:   1 = exclude the Condon-Shortley phase factor; -1 = include it.

lmax_calc : optional, integer, default = L
:   The maximum degree to recover. It must be between 0 and L.

# Description

`SHExpandCCC` independently projects positive and negative longitude orders
with parity-correct FFT interpolation and direct Clenshaw-Curtis quadrature.
It supports the maximum recoverable degree L=nlat-2, includes both poles, has
no `sampling` parameter, and does not use a latitude least-squares solve.
Unnormalized associated Legendre functions are accurate only to about degree
15.

# References

Clenshaw, C. W., and A. R. Curtis, *Numerische Mathematik*, 2, 197-205,
1960, https://doi.org/10.1007/BF01386223; Waldvogel, J., *BIT Numerical
Mathematics*, 46, 195-202, 2006, https://doi.org/10.1007/s10543-006-0045-4;
and Holmes, S. A., and W. E. Featherstone, *J. Geodesy*, 76, 279-299, 2002,
https://doi.org/10.1007/s00190-002-0216-2.
