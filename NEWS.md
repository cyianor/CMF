# CMF 1.0.4

## Internal changes
* Split indices from data values; this way indices can be stored as integers
  and data as numeric/floating point values. Theoretically reduces the storage
  need from n * 8 * 3 to n * (8 + 4 * 2), so 33% less. (WIP; needs to be
  implemented also from the start and not just internally of the estimation
  function)

## Minor changes
* Added a `NEWS.md` file to track changes to the package.
* Added `SystemRequirements: C++11` to `DESCRIPTION`
* Added `Language: en-US` to `DESCRIPTION`
* Removed `Date` from `DESCRIPTION`
* Add a `inst/CITATION` file for the original CMF article

# CMF 1.0.3

## Bug fixes
* Removed in-place manipulation of pseudo data, since that led to
  involuntary manipulation of the original data. Data is copied
  anyways, so creating a new matrix costs almost the same.

## Major changes
* Changed native code dependency from
  [Rcpp](https://cran.r-project.org/package=Rcpp)
  to [cpp11](https://cran.r-project.org/package=cpp1q)

## Minor changes
* Added a `README.md` file
* Added a `LICENSE.md` file
* Added `.Rbuildignore`
* Add bug tracker URL

## Miscellaneous
* Improved code formatting in R code
* Started using clang-format for C++ code
* Changed maintainer email

# CMF 1.0.2

## Bug fixes
* CRAN issues addressed to unarchive package

## Miscellaneous
* Change of maintainer
