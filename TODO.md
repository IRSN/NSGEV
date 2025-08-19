
### NSGEV TODO List

This is a temporary reminder list by YD. Not a planning for new
development!

- [ ] **`anova` name** Make an alias of the `anova` method under a
      better name making clear that this is for a Likelihood Ratio
      test. May be `LR.test`.
 
- [ ] **Profile Likelihood**. Speed up the computation of the profile
     likelihood intervals in `quantMax.TVGEV`. The initial values for
     the constrained optimisation could be better.

- [x] **Profile Likelihood**. Experiment on the computation of profile
     likelihood intervals with the ODE method? EXPERIMENTAL.
  
- [x] **Profile Likelihood**. Experiment on the computation of profile
      *contours* for some or all of the couples of parameters in a
      fitted `TVGEV` model.
  
- [] **Return level plots** Allow the return level plot to use an
      horizontal axis with a probability of exceedance (and a Gumbel
      scale) rather than a return period as done for the quantile of
      the maximum. So the plot for a conditional return level would be
      nearly identical to that given for the quantile of the maximum
      on a one-year period.
  
- [ ] **Improve Numerical Treatment**. Improve by using the QR method on
     the design matrices. So the vector of parameters will be
     internally $\mathbf{R}\boldsymbol{\psi}$ instead of
     $\boldsymbol{\psi}$. This should eliminate some of the problems
     met when the covariates are poorly scaled and/or are correlated.
  
- [x] **Tests**. Improve the code coverage by adding new tests. IN
      PROGRESS.
  
- [x] **Continuous Integration/Github Actions**. Add an action for the
     check, and an action for the build to produce a release when a
     tag is used on for the commit.
  
- [x] **Vignettes**. In the *Introduction* vignette enhance the part on
     the quantile of the maximum by adding profile likelihood
     intervals.  Complete the citations and introduce the expression
     *Design Life Level*.

- [x] **Vignettes**. In the *Introduction* vignette, include the fixes
     that were by Jesper Rydén (mail 2024-03-06).
  
- [ ] **README** file. Make the installation part shortest, and include
     a small example with plots that will be converted into a small
     Google Colab notebook.

- [ ] **Hexagon sticker**?

- [ ] Change the name of `bs` or replace this generic function by a
      generic bootstapping function from another package.

- [ ] Get rid of the dependence to the **reshape2** package by replacing
     `reshape2::melt` by `data.table::melt`.
