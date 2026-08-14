---
title: SHExpandCCC (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: shexpandccc.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Expand a complex Clenshaw-Curtis grid into spherical harmonics.

## Usage

call SHExpandCCC (`gridcc`, `cilm`, `lmax`, `norm`, `csphase`, `lmax_calc`,
`exitstatus`)

## Parameters

`gridcc` : input, complex(dp), dimension (`L`+2, 2\*`L`+2)
:   A non-extended complex Clenshaw-Curtis grid with equally spaced
    colatitudes and both poles.

`cilm` : output, complex(dp), dimension (2, `lmax_calc`+1, `lmax_calc`+1)
:   Complex spherical harmonic coefficients following the `SHExpandDHC`
    convention: positive orders in the first index and negative orders in the
    second before the conventional `(-1)^m` factor.

`lmax` : output, integer(int32)
:   The maximum recoverable degree inferred from the input: `L = nlat-2`.

`norm` : input, optional, integer(int32), default = 1
:   1 = geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized
    harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : input, optional, integer(int32), default = 1
:   1 = exclude the Condon-Shortley phase factor; -1 = include it.

`lmax_calc` : input, optional, integer(int32), default = `lmax`
:   The maximum degree to recover. It must be between 0 and `lmax`.

`exitstatus` : output, optional, integer(int32)
:   If present, return a status instead of stopping: 0 = no error; 1 = invalid
    dimensions; 2 = invalid input; 3 = allocation or FFTW-plan failure.

## Description

`SHExpandCCC` is the complex counterpart of `SHExpandCC`. It separately
interpolates and projects positive and negative longitude orders with the
same parity-correct FFT and direct Clenshaw-Curtis quadrature construction.
The input shape is (`L`+2, 2\*`L`+2), includes both poles, and supports the
maximum recoverable degree `L = nlat-2`. It has no `sampling` parameter and
does not perform a least-squares fit. Unnormalized associated Legendre
functions are accurate only to about degree 15.

## References

Clenshaw, C. W., and A. R. Curtis, A method for numerical integration on an
automatic computer, *Numerische Mathematik*, 2, 197-205, 1960,
https://doi.org/10.1007/BF01386223.

Waldvogel, J., Fast construction of the Fejer and Clenshaw-Curtis quadrature
rules, *BIT Numerical Mathematics*, 46, 195-202, 2006,
https://doi.org/10.1007/s10543-006-0045-4.

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
summation and the recursive computation of very high degree and order
normalised associated Legendre functions, *J. Geodesy*, 76, 279-299, 2002,
https://doi.org/10.1007/s00190-002-0216-2.

## See also

[makegridccc](makegridccc.html), [shexpandcc](shexpandcc.html),
[makegridcc](makegridcc.html), [shexpanddhc](shexpanddhc.html),
[shexpandglqc](shexpandglqc.html)
