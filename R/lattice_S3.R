# This file is part of the cry package
#

#########
### S3 classes related to the crystal lattice
########

#' Constructor for an S3 object of class "bravais"
#'
#' There are 14 Bravais lattices. They are represented by a two-letter symbol:
#' aP, mP, mS, oP, oS, oF, oI, tP, tI, hP, hR, cP, cF, cI
#'
#' @param bt A two-letter character, denoting the Bravais type.
#' @return An object of class "bravais". It is a named list of length 4. The first slot, "bt", is
#'         the universally-used two-letter symbol. The second, third and fourth slots are, respectively,
#'         "cr_fam" (the corresponding crystal family), "cr_sys" (the corresponding crystal system) and
#'         "lt_sys" (the corresponding lattice system).
#' @examples
#' # mS is a monoclinic, face-centred Bravais lattice
#' bt <- bravais("mS")
#' class(bt)
#' bt[1:4]
#'
#' @export
bravais <- function(bt=NULL) {
  # When no information is given the type is "aP"
  if (is.null(bt)) bt <- "aP"

  # Check input is a character
  if (!is.character(bt)) stop('Input is a two-letter character string.'
                              )
  # Check that input character has length 2
  if (nchar(bt) != 2) stop("The number of characters differs from 2")

  # Extract lattice and centering (also check on symbol validity)
  latt <- substr(bt,1,1)
  centring <- substr(bt,2,2)
  if (latt != "a" & latt !="m" & latt != "o" & latt != "t" & latt != "h" & latt != "c")
    stop("Invalid lattice type")
  if (centring != "P" & centring != "A" & centring != "B" & centring != "C" & centring != "F" &
      centring != "S" & centring != "I" & centring != "R") stop("Invalid centring type")

  # Now check validity for specific families
  # Triclinic family
  if (latt == "a" & centring != "P") stop("Invalid Bravais lattice type")

  # Monoclinic family
  if ((latt == "m" & centring == "F") | (latt == "m" & centring == "R"))
    stop("Invalid Bravais lattice type")

  # Orthorombic family
  if ((latt == "o" & centring != "P") & (latt == "o" & centring != "A") &
      (latt == "o" & centring != "B") & (latt == "o" & centring != "C") &
      (latt == "o" & centring != "S") &
      (latt == "o" & centring != "F") & (latt == "o" & centring != "I"))
    stop("Invalid Bravais lattice type")

  # Tetragonal family
  if ((latt == "t" & centring != "P") & (latt == "t" & centring != "I"))
    stop("Invalid Bravais lattice type")

  # Hexagonal and Trigonal families
  if (latt == "h"  & (centring != "P" & centring != "R")) stop("Invalid Bravais lattice type")

  # Cubic family
  if (latt == "c" & (centring == "A" | centring == "B" | centring == "C"))
    stop("Invalid Bravais lattice type")

  # All checks passed. Proceed to determine family, crystal system and lattice system
  if (latt == "a")
  {
    cr_fam <- "triclinic"
    cr_sys <- "triclinic"
    lt_sys <- "triclinic"
  }
  if (latt == "m")
  {
    cr_fam <- "monoclinic"
    cr_sys <- "monoclinic"
    lt_sys <- "monoclinic"
  }
  if (latt == "o")
  {
    cr_fam <- "orthorombic"
    cr_sys <- "orthorombic"
    lt_sys <- "orthorombic"
  }
  if (latt == "t")
  {
    cr_fam <- "tetragonal"
    cr_sys <- "tetragonal"
    lt_sys <- "tetragonal"
  }
  if (latt == "c")
  {
    cr_fam <- "cubic"
    cr_sys <- "cubic"
    lt_sys <- "cubic"
  }
  if (latt == "h")
  {
    cr_fam <- "hexagonal"
    if (centring == "R")
    {
      cr_sys <- "trigonal"
      lt_sys <- "rhombohedral or hexagonal (centred)"
    }
    if (centring == "P")
    {
      cr_sys <- "hexagonal"
      lt_sys <- "hexagonal"
    }
  }

  # S3 class
  bt <- structure(list(bt=bt,cr_fam=cr_fam,cr_sys=cr_sys,lt_sys=lt_sys),class = "bravais")

  return(bt)
}


#' Print method for an object of class "bravais".
#'
#' The Bravais lattice and related crystal family, crystal system and lattice system
#' are displayed.
#'
#' @param x An object of class "bravais".
#' @param ... Additional arguments passed to the print methods
#' @return Nothing. A message is displayed which includes
#'         information on the Bravais lattice.
#'
#' @examples
#' # Create a triclinic Bravais lattice
#' bt <- bravais()
#'
#' # Display its value
#' print(bt)
#' @rdname print.bravais
#' @export
print.bravais <- function(x,...) {
  cat("This is an object of class 'bravais'\n")
  msg <- sprintf("The Bravais lattive is:  %s.\n",x$bt)
  cat(msg)
  msg <- sprintf("Its crystal family is %s.\n",x$cr_fam)
  cat(msg)
  msg <- sprintf("Its crystal system is %s.\n",x$cr_sys)
  cat(msg)
  msg <- sprintf("Its lattice system is %s.\n",x$lt_sys)
  cat(msg)

  invisible(x)
}
