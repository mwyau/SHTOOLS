# MakeGridCC

Create a real Clenshaw-Curtis grid from spherical harmonic coefficients.

# Usage

call MakeGridCC (`gridcc`, `cilm`, `lmax`, `norm`, `csphase`, `lmax_calc`,
`extend`, `exitstatus`)

# Parameters

`gridcc` : output, real(dp), dimension (`lmax`+2, 2\*`lmax`+2+`extend`)
:   A Clenshaw-Curtis grid with both poles. When `extend` is 1, the last
    column is the redundant 360 E longitude.

`cilm` : input, real(dp), dimension (2, `lmaxin`+1, `lmaxin`+1)
:   Real spherical harmonic coefficients: cosine coefficients in the first
    index and sine coefficients in the second.

`lmax` : input, integer(int32)
:   The maximum degree defining the output grid.

`norm` : input, optional, integer(int32), default = 1
:   1 = geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized
    harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : input, optional, integer(int32), default = 1
:   1 = exclude the Condon-Shortley phase factor; -1 = include it.

`lmax_calc` : input, optional, integer(int32), default = `lmax`
:   The maximum degree synthesized. It must be between 0 and `lmax`.

`extend` : input, optional, integer(int32), default = 0
:   If 1, append the redundant longitude column at 360 E. No latitude row is
    appended because both poles are intrinsic to the grid.

`exitstatus` : output, optional, integer(int32)
:   If present, return a status instead of stopping: 0 = no error; 1 = invalid
    dimensions; 2 = invalid input; 3 = allocation or FFTW-plan failure.

# Description

`MakeGridCC` evaluates associated Legendre functions at `lmax+2` equally
spaced colatitudes from 0 through pi and uses an FFT over `2*lmax+2`
longitudes. It is the inverse of `SHExpandCC`. The non-extended grid has
shape (`lmax`+2, 2\*`lmax`+2), includes both poles, and has no `sampling`
parameter. The Legendre recurrences follow Holmes and Featherstone (2002);
unnormalized functions are accurate only to about degree 15.

# Reference

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
summation and the recursive computation of very high degree and order
normalised associated Legendre functions, *J. Geodesy*, 76, 279-299, 2002,
https://doi.org/10.1007/s00190-002-0216-2.

# See also

[shexpandcc](shexpandcc.html), [makegridccc](makegridccc.html),
[shexpandccc](shexpandccc.html), [makegriddh](makegriddh.html),
[makegridglq](makegridglq.html)
