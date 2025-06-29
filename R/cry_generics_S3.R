# This file is part of the cry package
#

#########
### S3 generics specific for cry objects
#
# !!! All @imports in create_unit_cell are not just for unit cell, but to make Roxygen2
# add them automatically to DESCRIPTION (basically I had to add them somewhere!)
########

#' S3 generic to create unit_cell objects
#'
#' The unit_cell object can be created starting from specific objects, files, etc.
#'
#' @param a An object or objects used to select a method. These can be
#'          cell parameters, an object of class rec_unit_cell, etc.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class "unit_cell". It is a named list of length 6 whose
#'         last three slots are of class "angle".
#'
#' @examples
#' # Create a unit_cell in default (no arguments)
#' uc <- create_unit_cell()
#' print(uc)
#' @importFrom graphics hist
#' @importFrom stats na.omit sd
#' @importFrom utils read.table write.table
#' @importFrom methods is
#' @importFrom stats complete.cases fft
#'
#' @export
create_unit_cell <- function(a,...) UseMethod("create_unit_cell")


#' S3 generic to create rec_unit_cell objects
#'
#' The rec_unit_cell object can be created starting from specific objects, files, etc.
#'
#' @param ar An object or objects used to select a method. These can be
#'           reciprocal unit cell parameters, an object of class
#'           rec_unit_cell, etc.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class "rec_unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#'
#' @examples
#' # Create a rec_unit_cell in default (no arguments)
#' ruc <- create_rec_unit_cell()
#' print(ruc)
#' @export
create_rec_unit_cell <- function(ar,...) UseMethod("create_rec_unit_cell")


#' S3 generic to compute cell volume
#'
#' The volume of a unit cell and a reciprocal unit cell can be calculated starting
#' from specific objects, files, etc.
#'
#' @param x An object used to select a method. Either an object of class
#'          unit_cell or an object of class rec_unit_cell.
#' @param ... Further arguments passed to or from other methods.
#' @return V A real number. The volume (in \eqn{angstroms^3} or
#'           \eqn{angstroms^{-3}}) of the input unit cell or reciprocal
#'           unit cell.
#'
#' @examples
#' # Calculate the volume of a unit cell
#' uc <- unit_cell(20)
#' V <- calculate_cell_volume(uc)
#'
#' # Calculate the volume of the corresponding reciprocal cell
#' ruc <- create_rec_unit_cell(uc)
#' Vrec <- calculate_cell_volume(ruc)
#' V*Vrec  # Should be 1!
#'
#' @export
calculate_cell_volume <- function(x,...) UseMethod("calculate_cell_volume")

#' S3 generic to create merged_reflections objects
#'
#' The merged_reflections object can be created starting from
#' specific objects, files, etc.
#'
#' @param ruc An object used to select a method.
#' @param ... Further arguments passed to or from other methods.
#' @return mrefs An object of class "merged_reflections". It is a named
#'         list of length 4 whose names are:
#'         \describe{
#'           \item{ruc}{An object of class "rec_unit_cell".}
#'           \item{csym}{An object of class "cryst_symm".}
#'           \item{records}{A data frame containing the data.}
#'           \item{dtypes}{A character vector containing the
#'                 type of data (Miller indices, structure
#'                 factors, etc).}
#'         }
#'
#' @examples
#' # Create a default merged_reflections object (no arguments)
#' mrefs <- create_merged_reflections()
#' print(mrefs)
#'
#' # Create merged_reflections object from symmetry
#' csym <- cryst_symm("P 3")
#' mrefs <- create_merged_reflections(csym=csym)
#' print(mrefs)
#'
#' @importFrom graphics hist
#' @importFrom stats na.omit sd
#' @importFrom utils read.table write.table
#' @importFrom methods is
#' @export
create_merged_reflections <- function(ruc,...)
  UseMethod("create_merged_reflections")
