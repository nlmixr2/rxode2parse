# Cran comments

- This fixes some bugs and is required for the upcoming StanHeaders release on CRAN

- This require C++14 because of StanHeaders

- This also updates the URLs of codecov.io

- Adding this new version of `rxode2parse` will break `rxode2` on
  CRAN.  The following packages will be updated to support a new
  version of `rxode2`: `rxode2et` and `rxode2random`.

- `nlmixr2est` has a binary linkage to `rxode2` and will also be
  updated.
