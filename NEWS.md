# HonestDiD 0.2.9

- Fix remaining `CVXR::psolve()` calls in `.maxBiasFN` and `.minBiasFN` in `arp-nuisance.R` to use `.psolve()` wrapper for CVXR >= 1.8 compatibility

# HonestDiD 0.2.8

- Restore legacy CVXR::psolve interface (`$solve`, `$status`, `$getValue`)
- Add `tests/test_syntax.R` with better options coverage

# HonestDiD 0.2.7

* Replace `CVXR::solve()` with `CVXR::psolve()` for compatibility with CVXR >= 1.0,
  which no longer exports `solve` (#68).

# HonestDiD 0.2.6

* Improve error messages; implemented scaling retry for `.findLowestH` (#63).
