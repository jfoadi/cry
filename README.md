
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cry <img src='man/figures/logo.png' align="right" height="128" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/cry)](https://CRAN.R-project.org/package=cry)
<!-- badges: end -->

*cry* is an R package to make it easier dealing with crystallographic
data. The package includes functions to read/write data from/to some of
the file mostly used in software for structural crystallography. Current
entry include the `readMTZ` function to read MTZ files (see
[CCP4](https://www.ccp4.ac.uk)) and the corresponding `writeMTZ`. There
are also `readCIF`, `readpd_rtv` and `readSF_CIF` to read files in CIF
formats (Crystallographic Information Framework, see [IUCr
page](https://www.iucr.org/resources/cif)). Users are welcome to suggest
inclusions of different (not currently available) formats they might
need for specific tasks and/or analysis.

*cry* includes also several functions to perform the most common and
routine crystallographic calculations, many of them involving
crystallographic symmetry. The main purpose of these functions is to
enable users to investigate specific issues without having to resort to
external packages and thus carrying the analysis through without having
to abandon the R platform.

## Installation

You can install the released version of crone from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("cry")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jfoadi/cry")
```

## Example

This is a simple example in which a CIF file containing structure
factors is read and a part of its data used for statistical analysis.
The file is included as example file in *cry* and is called
*1dei-sf.cif*. The first task is to import the data in R.

``` r
# Load cry package
library(cry)

# Save CIF data in a named list
datadir <- system.file("extdata",package="cry")
filename <- file.path(datadir,"1dei-sf.cif")
lCIF <- readSF_CIF(filename)
```

The newly-created R object, a list called `lCIF`, contains all that is
included in the CIF file, but within a structure that can be used within
R for various analyses. The components of this list can be readily
explored using the common R functions `names`, `class` and `str`.

``` r
# What's containe in lCIF?
names(lCIF)
#> [1] "HEADER" "SYMM"   "REFL"

# What's the header?
class(lCIF$HEADER)   # It's a list
#> [1] "list"

# Is it a named list?
names(lCIF$HEADER)  # Yes
#> [1] "TITLE"      "CELL"       "SGN"        "HALL"       "HM"        
#> [6] "SHELX_CODE" "F_000"      "HIGH_RES"

# What's the space group name for this crystal structure?
print(lCIF$HEADER$HM)
#> [1] NA

# What's the space group number corresponding to P 21 21 21?
# Use one of cry's functions
xHM <- lCIF$HEADER$HM
translate_SG(xHM,SG_in="xHM",SG_out="number")$msg
#> [1] "Something wrong in your input:\n   1) the symbol or number input for this space group does not exist\n   2) if your inpur was a number, perhaps for this space group there are not that many settings"

# ... and the unit cell parameters?
cpars <- c(lCIF$HEADER$CELL$A$VAL,lCIF$HEADER$CELL$B$VAL,
           lCIF$HEADER$CELL$C$VAL,lCIF$HEADER$CELL$ALPHA$VAL,
           lCIF$HEADER$CELL$BETA$VAL,lCIF$HEADER$CELL$GAMMA$VAL)
print(cpars)
#> [1] NA NA NA NA NA NA

# The unit cell belongs to the orthorombic system,
# as it should be, due to symmetry (cry function)
print(crystal_system(gn=19))
#> [1] "ORTHOROMBIC"
```

The most important bit in the CIF file are the observed reflections
(coming from x-ray diffraction of one or more crystals). They are
contained in the named list `REFL`. All CIF files have a `VAL` field
reporting the specific observed quantity and a field `STD` reporting the
corresponding experimental error. Not always the last one is available.
The `VAL` field is an R data frame that makes these data suitable to
further analysis within R. The `STD` field mirrors `VAL`, but with the
experimental errors counterparts, if available, otherwise it is `NULL`,
as in the case here reported.
