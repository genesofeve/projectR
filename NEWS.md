# projectR 1.24.0

## Bug fixes

- `projectR()` prcomp, matrix, dgCMatrix, and rotatoR dispatches no longer
  drop the matrix dimension when `NP` is a scalar integer. Previously, a single
  integer `NP` caused `loadings[, NP]` to silently return a vector, producing
  an `incorrect number of dimensions` error inside `geneMatchR` (#50).

## New features

- `projectR()` matrix, dgCMatrix, and LinearEmbeddingMatrix (CoGAPS/NMF)
  dispatches gain an `include_intercept` argument (default `FALSE`). When
  `TRUE`, the design is fit as `~ 1 + A` instead of `~ 0 + A`, adding an
  explicit intercept that absorbs the per-sample global mean shift and leaves
  pattern coefficients to describe only pattern-specific variation. When
  `full = TRUE`, the per-sample intercept estimates are returned in
  `$intercept` (#51).

- `projectR()` prcomp dispatch gains a `center_by_loadings` argument. When
  `TRUE`, the new data is centered using the per-gene means stored in
  `loadings$center` (computed from the original training data) rather than by
  its own column means. This is useful when projecting longitudinal or paired
  samples where mean-centering by the new data would remove biologically
  meaningful shifts (#46).

- `projectR()` `full = TRUE` now returns a per-sample coefficient of
  determination (`$r_squared`) for both the matrix/NMF/CoGAPS and prcomp
  dispatch methods. R-squared is computed as
  `1 - ||residual||^2 / ||X||^2` using the raw (uncentered) total sum of
  squares, which keeps values comparable across dispatch methods (#43).
  
- `NEWS.md` added to track changes/updates

## Improvements

- `projectR()` matrix dispatch p-values now use the t-distribution with
  `lmFit` residual degrees of freedom instead of the standard normal
  approximation. The previous approximation was anti-conservative when the
  number of samples was small relative to the number of patterns; the
  t-distribution is exact for all sample sizes and converges to the normal
  for large degrees of freedom (#49).

- `projectR()` prcomp dispatch `full = TRUE` now returns a named list with
  elements `$projection`, `$pvar`, and `$r_squared`, consistent with the
  naming convention used by the matrix/NMF/CoGAPS dispatch. Previously the
  prcomp method returned an unnamed two-element list (#48).

- `?projectR` Docstring Value section now uses subsection headings to
  distinguish the prcomp dispatch return value (which contains `$pvar`)
  from the matrix/NMF/CoGAPS dispatch return value (which contains `$pval`).

- `multivariateAnalysisR()` updated to replace the deprecated `geom_errorbarh()`
  with `geom_errorbar(orientation = "y")` and `size` aesthetic with `linewidth`,
  resolving ggplot2 (>= 3.4.0 / 4.0.0) deprecation warnings.

- `R/projectR.R` styled with `styler` for a more consistent code appearance.

- `R/.cluster2pattern.R` was documented with roxygen comments, but not exported.
  Raised a roxygen warning; documentation changed to `#` to avoid the issue.

## Tests

- Added tests for `projectR()` prcomp dispatch `full = TRUE`: verifies the
  returned list has names `$projection`, `$pvar`, and `$r_squared`, and that
  `$r_squared` values are in range (#43, #48).

- Added tests for `projectR()` matrix dispatch `full = TRUE`: verifies the
  returned list has names `$projection`, `$pval`, and `$r_squared`, and that
  `$pval` values are in [0, 1] (#43, #49).

- Added tests for `center_by_loadings`: verifies the two centering paths
  produce different projections, and that passing `center_by_loadings = TRUE`
  with a prcomp fit without centering raises an informative error (#46).

- Added test for `projectR()` prcomp dispatch with scalar `NP`: verifies that
  a single integer `NP` with `center_by_loadings = TRUE` and `full = TRUE`
  does not error and returns a projection with one row (#50).

- `multivariateAnalysisR()` test now uses `withr::local_dir(withr::local_tempdir())`
  to redirect file output to a temporary directory, preventing CSV and PNG
  files from being written into the test directory. File existence is verified
  explicitly before the temp directory is cleaned up.

# projectR 1.23.2

- `projectR()` sparse matrix (dgCMatrix) dispatch `full = TRUE` now returns a
  named list consistent with the dense matrix dispatch format (#44).

# projectR 1.23.1

- `projectR()` gains support for sparse matrix (dgCMatrix) input (#39). A
  `chopBy` argument controls chunked projection for large datasets.

- Fix missing class definition warning for `LinearEmbeddingMatrix`.

- Fix broken Bioconductor vignette link (#41).
