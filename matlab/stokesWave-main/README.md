# stokesWave

Reference article: Zhao and Liu (2022), On Stokes Wave Solutions. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, Volume 478, Issue 2258, https://doi.org/10.1098/rspa.2021.0732. 

List of codes.

1. StokesTheoryCombined.mlx 

This file is used to derive solutions for Stokes waves following the manuscript 'On Stokes Wave Solutions' by Zhao and Liu (2022). To be used with Matlab Live Editor for best view of equations. Note that running this file consumes lot of computer memory, and the fifth order derivation takes relatively long time. Try to close other programmes and save your work if you intend to run this file in live editor. In this file, the symbol H means the wave height for first order first harmonic term for derivation. This is not the same H as in the manuscript. 

2. toUseStokesTheory.m 

This file demonstrate how to use the functions prepared for Stokes wave.

3. StokesDispSolver.m 

This function is a solver for nonlinear dispersion relationship. Different cases were provided to find wave length, for example (1) given a, h, T, find H and k; (2) given H, h, T, find a and k; (3) given aw, h, T, find Hw and k, etc.

4. StokesEta.m 

With the wave amplitude and wave number solved using Function 3, calculate the free surface using Stokes theory.

5. StokesU.m

With the wave amplitude and wave number solved using Function 3, calculate the velocity using Stokes theory.

6. normalizeFreeSurface.mlx

This file is to be opened using Matlab Live Editor. This file shows the process of converting the free surface equations in S and H, and Fenton to the format in the present manuscript for comparison. 

7. normalizePotentialFunction.mlx

This file shows the process of converting the potential functions in S and H, and Fenton to the format in the present manuscript for comparison. 

8. OtherSolutions.zip

This group of solution is the function to find results using Fenton's equations.

FentonU.m

FentonEta.m

FentonDispSolver.m

This group of solution is the function to find results using Skjelbreia and Hendrickson's equations.

SnHDispSolver.m

SnHEta.m

SnHU.m

Others Codes

hyperIden.m

This is a function to change format of hyperbolic function to \cosh () using the hyperbolic identities.

trigIden.m

This is a function to change format of trigonometric function to \sin () or \cos () using the trigonometric identities.

