# SHExpandCC

Expand a real Clenshaw-Curtis grid into spherical harmonics.

# Usage

call SHExpandCC (`gridcc`, `cilm`, `lmax`, `norm`, `csphase`, `lmax_calc`,
`exitstatus`)

# Parameters

`gridcc` : input, real(dp), dimension (`L`+2, 2\*`L`+2)
:   A non-extended Clenshaw-Curtis grid with equally spaced colatitudes and
    both poles.

`cilm` : output, real(dp), dimension (2, `lmax_calc`+1, `lmax_calc`+1)
:   The real cosine and sine spherical harmonic coefficients.

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

# Description

`SHExpandCC` performs a direct regular-grid analysis: FFTs transform longitude
orders, parity-correct continuation and zero-padding interpolate each order to
a twice-finer colatitude grid, and classical Clenshaw-Curtis quadrature
projects onto associated Legendre functions. The fine grid has 2\*`L`+3
latitude nodes, which integrate the degree-2\*`L` products in the projection.
No latitude least-squares solve or dense latitude-by-degree matrix is built.

The input shape is (`L`+2, 2\*`L`+2); both poles are retained and the maximum
recoverable degree is `L`. CC has one layout and no `sampling` parameter. The
normalized Legendre recurrences follow Holmes and Featherstone (2002);
unnormalized functions are accurate only to about degree 15.

# References

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

# See also

[makegridcc](makegridcc.html), [shexpandccc](shexpandccc.html),
[makegridccc](makegridccc.html), [shexpanddh](shexpanddh.html),
[shexpandglq](shexpandglq.html)
