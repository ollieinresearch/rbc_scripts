# Changing resolution exactly

{mod}`rbc.regrid` changes the resolution of a coefficient-space checkpoint by padding
every axis with zeros (to refine) or truncating it (to coarsen). Refining is exact; coarsening
removes exactly the modes the coarser grid cannot hold.

## Why zero-padding is exact

A field is stored as coefficients of fixed basis functions:

$$
f(x, y, z) = \sum_{k_x, k_y, n} c_{k_x k_y n}\;\phi_{k_x}(x)\,\phi_{k_y}(y)\,p_n(z),
$$

with $\phi_k$ the unit-amplitude cosines and sines of mode $k$ and $p_n$ the orthonormal
Chebyshev polynomials. The key fact, checked in the Dedalus source
({mod}`rbc.spectra.core`), is that **neither the basis functions nor the normalization of
the coefficients depends on the number of modes**: Dedalus divides by $N$ in its forward
transforms. A field with $N_z$ modes is therefore the same function as the field with
$2N_z$ modes whose first $N_z$ coefficients are the same and the rest zero.

Coarsening keeps the first coefficients and drops the rest, which is the best
approximation of the field in the smaller space (orthogonal projection in the
coefficients' inner product). The dropped part is what a coarser run cannot represent
anyway, and its size is visible beforehand in the spectra.

## What else is handled

- The Fourier storage interleaves cosine and sine coefficients, two entries per mode
  ($k$ at entries $2k$ and $2k+1$), and has no Nyquist entry. Keeping the first $N_x$
  entries keeps modes $0 \dots N_x/2 - 1$, so $N_x$ must be even.
- Tau fields depend only on $(x, y)$: their size-1 $z$ axis is left alone.
- The output holds one write, as plain datasets, which Dedalus' own loader reads on any
  process mesh.
- It is written to a temporary name and renamed at the end, so a failure never leaves a
  half-written restart file.

A grid-space checkpoint needs none of this: Dedalus interpolates those itself when the
shape changes, and {mod}`rbc.regrid` refuses them.

## How it was checked

The implementation was tested with Dedalus itself: a run with $N_z = 16$ was regridded to
$N_z = 32$ and continued, and the volume averages just after the restart matched those of
the coarse run to a relative difference of about $10^{-6}$.
