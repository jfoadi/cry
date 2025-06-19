## Submission of version 0.5.2
Triggered by ggplot2 being changed. Also, sorted a few small issues.
## Submission of version 0.5.1
This is the first release after a long time. Prompted mostly by request due to a failure from CRAN and M. Maechler (for R-devel).
With this version I have also experienced GitHub Actions for the first time.

## Submission of version 0.5.0
After first release in CRAN, this is the version containing many more functions and improvements.


## Resubmission
This is a resubmission. Compared to the first submission attempt, I have:

* Added a reference to DESCRIPTION and updated the date.

* Included \value to some methods I previously missed.

* Uncommented and changed some examples who were previously commented out.

* Eliminated a couple of methods as I deemed them incomplete.


## Test environments
* Windows 10 - Using RStudio 1.2.5033
* Linux - Using TRAVIS

## R CMD check results
Carried out from within RStudio.
No ERRORS.
1 WARNING:
 WARNING
'qpdf' is needed for checks on size reduction of PDFs
This doesn't pop up in the building carried out by Travis and with devtools::check_rhub.
No NOTES.
This is the first resubmission, after initial submission.

## Downstream dependencies
There are currently no downstream dependencies for this package.
