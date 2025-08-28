# Introduction

Provided are MATLAB source codes to compute the Molodensky's truncation
coefficients from the [Bucha et
al. (2019a)](https://doi.org/10.1007/s00190-018-1139-x) and [Bucha et
al. (2019b)](https://doi.org/10.1007/s00190-019-01277-3) studies.  The
coefficients can be used to limit the integration domain in spectral gravity
forward modelling to near- and/or far-zone masses.  The
[ADVANPIX](https://www.advanpix.com/) toolbox is required to run the code.

For a test run, see the documentation in the source code files inside the `src`
folder.


## Pre-computed truncation coefficients

Some pre-computed ready-to-use near-zone and far-zone truncation coefficients
from the [Bucha et al. (2019b)](https://doi.org/10.1007/s00190-019-01277-3)
study are available at
[https://doi.org/10.5281/zenodo.7092047](https://doi.org/10.5281/zenodo.7092047).


# Important note

Prior to the `de6f4916f967d39edbd092274e35b5a2bcb6ff34` commit, the Q11, Q21
and Q22 coefficients were computed *incorrectly* for topography powers larger
than 3.  Do not use the old code.

When downloading the pre-computed coefficients from the Zenodo link above,
always make sure that you are downloading version 2 or newer.  Do not use the
coefficients from version 1, as they are wrong for Q11, Q21 and Q22.


# Recommendation

This code and this repository are deprecated.  Use
[CHarm](https://www.charmlib.org) instead.  It is orders of magnitude faster
and relies solely on free software.  CHarm can be called from C and Python.


# Contact

Feel free to contact the author, Blazej Bucha, at blazej.bucha@stuba.sk.


# Citing

* Bucha, B., Hirt, C., Kuhn, M., 2019. Cap integration in spectral gravity
  forward modelling: near- and far-zone gravity effects via Molodensky's
  truncation coefficients. Journal of Geodesy 93, 65-83,
  [https://doi.org/10.1007/s00190-018-1139-x](https://doi.org/10.1007/s00190-018-1139-x)

* Bucha, B., Hirt, C., Kuhn, M., 2019. Cap integration in spectral gravity
  forward modelling up to the full gravity tensor. Journal of Geodesy 93,
  1707-1737,
  [https://doi.org/10.1007/s00190-019-01277-3](https://doi.org/10.1007/s00190-019-01277-3)


# Other related projects

* [CHarm](https://github.com/blazej-bucha/charm): C library for spherical
  harmonic transforms up to high degrees (tens of thousands and beyond).
  Supports OpenMP parallelization for shared memory architectures and
  vectorized CPU instructions (AVX, AVX2, AVX-512).

* [GrafLab](https://github.com/blazej-bucha/graflab) (GRAvity Field
  LABoratory): GrafLab (GRAvity Field LABoratory) is a MATLAB-based routine to
  compute functionals of the geopotential up to high degrees (tens of thousands
  and beyond).

* [isGrafLab](https://github.com/blazej-bucha/isgraflab) (Irregular Surface
  GRAvity Field LABoratory): A modified version of GrafLab to perform accurate
  and fast synthesis of gravity field quantities at dense grids residing on
  irregular surfaces such as the Earth's surface.

