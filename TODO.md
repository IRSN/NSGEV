
### NSGEV TODO List

This is a temporary reminder list by YD. Not a planning for new
development!

- [] **Profile Likelihood**. Speed up the computation of the profile
  likelihood intervals in `quantMax.TVGEV`. The initial values for the
  constrained optimisation could be better.

- [] **Profile Likelihood**. Experiment on the computation of profile
  likelihood intervals with the ODE method?
  
- [] **Improve Numerical Treatment**. Improve by using the QR method
  on the design matrices. So the vector of parameters will be
  internally $\mathbf{R}\boldsymbol{\psi}$ instead of
  $\boldsymbol{\psi}$. This should eliminate some of the problems met
  when the covariates are poorly scaled and/or are correlated.
  
- [] **Tests**. Improve the code coverage by adding new tests.
  
- [] **Continuous Integration/Github Actions**. Add an action for the
  check, and an action for the build to produce a release when a tag is
  used on for the commit.
  
- [] **Vignettes**. In the Introduction vignette enhance the part on
  the quantile of the maximum by adding profile likelihood intervals.
  Complete the citations and introduce the expression *Design Life
  Level*.

- [] **Vignettes**. In the Introduction vignette, include the fixes
  that were by Jesper Rydén (mail 2024-03-06).
  
- [] **README** file. Make the installation part shortest, and include
  a small example with plots that will be converted into a small
  Google Colab notebook.

- [] **Hexagon sticker**?
