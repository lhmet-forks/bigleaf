## Test environments
- x86_64-w64-mingw32/x64 (64-bit); R version 3.5.1 (package development)
- checks on other environments using the rhub package:
  - Debian Linux, R-release, GCC
  - Debian Linux, R-devel, GCC ASAN/UBSAN
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results
0 errors, 0 warnings, 0 notes

## Downstream dependencies
the package has the following revers dependencies:
Reverse suggests: REddyProc

## Resubmission
this is a resubmission in response to the email by Kurt Hornik from 26/05/2019.
The package was rebuilt with an updated version (1.23) of the knitr package.
