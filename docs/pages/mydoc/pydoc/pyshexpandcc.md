---
title: SHExpandCC()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshexpandcc.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Expand a real Clenshaw-Curtis grid into spherical harmonics.

## Usage

cilm = SHExpandCC(gridcc, [norm, csphase, lmax_calc])

## Returns

cilm : float, dimension (2, lmax_calc+1, lmax_calc+1)
:   Real cosine and sine spherical harmonic coefficients.

## Parameters

gridcc : float, dimension (L+2, 2\*L+2)
:   A non-extended Clenshaw-Curtis grid containing both poles.

norm : optional, integer, default = 1
:   1 = geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized
    harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

csphase : optional, integer, default = 1
:   1 = exclude the Condon-Shortley phase factor; -1 = include it.

lmax_calc : optional, integer, default = L
:   The maximum degree to recover. It must be between 0 and L.

## Description

`SHExpandCC` analyzes through the maximum recoverable degree L=nlat-2 without
discarding either pole. It uses FFTs in longitude, parity-correct spectral
interpolation in colatitude, and direct Clenshaw-Curtis quadrature projection
onto the associated Legendre basis. It is not a Driscoll and Healy transform
and does not use a latitude least-squares solve. CC has no `sampling`
parameter. Unnormalized associated Legendre functions are accurate only to
about degree 15.

## References

Clenshaw, C. W., and A. R. Curtis, *Numerische Mathematik*, 2, 197-205,
1960, https://doi.org/10.1007/BF01386223; Waldvogel, J., *BIT Numerical
Mathematics*, 46, 195-202, 2006, https://doi.org/10.1007/s10543-006-0045-4;
and Holmes, S. A., and W. E. Featherstone, *J. Geodesy*, 76, 279-299, 2002,
https://doi.org/10.1007/s00190-002-0216-2.
