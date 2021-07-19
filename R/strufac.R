#
# This file is part of the cry package
#

#' Coefficients for scattering factors
#'
#' The nine coefficients of the four-gaussians approximation for
#' the scattering factors
#'
#' A scattering factor is a function of
#' \eqn{\sigma\equiv\sin(\theta)/\lambda},
#' which is equal to \eqn{s/2}, with \eqn{s=1/d} and \eqn{d} the
#' resolution of the given reflection. The scattering form is
#' satisfactorily modelled as sum of four gaussians:
#' \deqn{
#' a_1\exp(-b_1 \sigma^2) + a_2\exp(-b_2 \sigma^2) +
#' a_3\exp(-b_3 \sigma^2) + a_4\exp(-b_4 \sigma^2) + c
#' }
#' This function returns the coefficients as a named list with
#' three components, \eqn{a=(a_1,a_2,a_3,a_4)},
#' \eqn{b=(b_1,b_2,b_3,b_4)} and \eqn{c=c}. The possible atom names
#' (also included are some ionic species) are here summarized.
#' \itemize{
#' \item H
#' \item H-1
#' \item He
#' \item Li
#' \item Li+1
#' \item Be
#' \item Be+2
#' \item B
#' \item C
#' \item Cv
#' \item N
#' \item O
#' \item O-1
#' \item F
#' \item F-1
#' \item Ne
#' \item Na
#' \item Na+1
#' \item Mg
#' \item Mg+2
#' \item Al
#' \item Al+3
#' \item Si
#' \item Siv
#' \item Si+4
#' \item P
#' \item S
#' \item Cl
#' \item Cl-1
#' \item Ar
#' \item K
#' \item K+1
#' \item Ca
#' \item Ca+2
#' \item Sc
#' \item Sc+3
#' \item Ti
#' \item Ti+2
#' \item Ti+3
#' \item Ti+4
#' \item V
#' \item V+2
#' \item V+3
#' \item V+5
#' \item Cr
#' \item Cr+2
#' \item Cr+3
#' \item Mn
#' \item Mn+2
#' \item Mn+3
#' \item Mn+4
#' \item Fe
#' \item Fe+2
#' \item Fe+3
#' \item Co
#' \item Co+2
#' \item Co+3
#' \item Ni
#' \item Ni+2
#' \item Ni+3
#' \item Cu
#' \item Cu+1
#' \item Cu+2
#' \item Zn
#' \item Zn+2
#' \item Ga
#' \item Ga+3
#' \item Ge
#' \item Ge+4
#' \item As
#' \item Se
#' \item Br
#' \item Br-1
#' \item Kr
#' \item Rb
#' \item Rb+1
#' \item Sr
#' \item Sr+2
#' \item Y
#' \item Y+3
#' \item Zr
#' \item Zr+4
#' \item Nb
#' \item Nb+3
#' \item Nb+5
#' \item Mo
#' \item Mo+3
#' \item Mo+5
#' \item Mo+6
#' \item Tc
#' \item Ru
#' \item Ru+3
#' \item Ru+4
#' \item Rh
#' \item Rh+3
#' \item Rh+4
#' \item Pd
#' \item Pd+2
#' \item Pd+4
#' \item Ag
#' \item Ag+1
#' \item Ag+2
#' \item Cd
#' \item Cd+2
#' \item In
#' \item In+3
#' \item Sn
#' \item Sn+2
#' \item Sn+4
#' \item Sb
#' \item Sb+3
#' \item Sb+5
#' \item Te
#' \item I
#' \item I-1
#' \item Xe
#' \item Cs
#' \item Cs+1
#' \item Ba
#' \item Ba+2
#' \item La
#' \item La+3
#' \item Ce
#' \item Ce+3
#' \item Ce+4
#' \item Pr
#' \item Pr+3
#' \item Pr+4
#' \item Nd
#' \item Nd+3
#' \item Pm
#' \item Pm+3
#' \item Sm
#' \item Sm+3
#' \item Eu
#' \item Eu+2
#' \item Eu+3
#' \item Gd
#' \item Gd+3
#' \item Tb
#' \item Tb+3
#' \item Dy
#' \item Dy+3
#' \item Ho
#' \item Ho+3
#' \item Er
#' \item Er+3
#' \item Tm
#' \item Tm+3
#' \item Yb
#' \item Yb+2
#' \item yb+3
#' \item Lu
#' \item Lu+3
#' \item Hf
#' \item Hf+4
#' \item Ta
#' \item Ta+5
#' \item W
#' \item W+6
#' \item Re
#' \item Os
#' \item Os+4
#' \item Ir
#' \item Ir+3
#' \item Ir+4
#' \item Pt
#' \item Pt+2
#' \item Pt+4
#' \item Au
#' \item Au+1
#' \item Au+3
#' \item Hg
#' \item Hg+1
#' \item hg+2
#' \item Tl
#' \item Tl+1
#' \item Tl+3
#' \item Pb
#' \item Pb+2
#' \item Pb+4
#' \item Bi
#' \item Bi+3
#' \item Bi+5
#' \item Po
#' \item At
#' \item Rn
#' \item Fr
#' \item Ra
#' \item Ra+2
#' \item Ac
#' \item Ac+3
#' \item Th
#' \item Th+4
#' \item Pa
#' \item U
#' \item U+3
#' \item U+4
#' \item U+6
#' \item Np
#' \item Np+3
#' \item Np+4
#' \item Np+6
#' \item Pu
#' \item Pu+3
#' \item Pu+4
#' \item Pu+6
#' \item Am
#' \item Cm
#' \item Bk
#' \item Cf
#' }
#'
#' @param atom_name A character string. The symbol of the atom whose
#'        scattering factor's coefficients are needed (see list of
#'        available names in details).
#' @return A named list whose three components are vectors and
#'         scalars containing the coefficients:
#'         \describe{
#'           \item{a}{Vector of length 4 containing the \eqn{a_j}}
#'           \item{b}{Vector of length 4 containing the \eqn{b_j}}
#'           \item{c}{Single number, the value of \eqn{c}}
#'         }
#'
#' @examples
#'
#' # Coefficients for the carbon atom
#' aname <- "C"
#' lcoeffs <- sfcoeffs(aname)
#' print(lcoeffs)
#'
#' @export
sfcoeffs <- function(atom_name) {
  # Make sure input name has correct length, or discard it
  atom_name <- trimws(atom_name,which="both")

  # Available names
  idx <- which(nchar(sfinfo) == 6)
  avnames <- trimws(sfinfo[idx],which="both")
  ans <- atom_name %in% avnames
  if (!ans) {
    msg <- "Input name is not in list available.\n"
    cat(msg)
    return(NULL)
  }

  # Format name into a 6-character string
  atom_name <- sprintf("%-6s",atom_name)

  # Extract interested-lines numbers from sfinfo
  istart <- which(sfinfo == atom_name)
  jc <- istart+1
  ja <- istart+2
  jb <- istart+3

  # A coefficients
  stmp <- strsplit(sfinfo[ja],"      ")[[1]]
  A <- na.omit(as.numeric(stmp))
  attributes(A) <- NULL

  # B coefficients
  stmp <- strsplit(sfinfo[jb],"      ")[[1]]
  B <- na.omit(as.numeric(stmp))
  attributes(B) <- NULL

  # C coefficient
  stmp <- strsplit(sfinfo[jc],"      ")[[1]]
  C <- as.numeric(stmp[length(stmp)])

  # Final named list
  lcoeffs <- list(a=A,b=B,c=C)

  return(lcoeffs)
}


#' Atomic scattering factor
#'
#' The atomic scattering factor function used in x-ray
#' crystallography
#'
#' The scattering factor is normally built as a function of
#' \eqn{\sigma\equiv\sin\theta/\lambda=s/2}, where \eqn{s}
#' is the length of the Miller indices (as a vector) in reciprocal
#' space, \eqn{s=|\mathbf{h}|}. But in the expression for the
#' structure factor it appears as a function of \eqn{s}, rather
#' than a function of \eqn{\sigma}. Here the scattering function
#' is a function of \eqn{s}. The approximation used here (and in
#' most crystallography packages) is accurate for up to 2 angstroms
#' resolution. It can be used for resolutions higher than 2
#' angstroms, but it is less accurate.
#'
#' @param atom_name A character string. The symbol of the atom whose
#'        scattering factor needs to be calculated (see list of
#'        available names in \code{\link{sfcoeffs}}).
#' @param s A real numeric vector, the values of \eqn{s} at which
#'           the scattering function needs to be evaluated.
#' @return The value of the scattering factor at s.
#'
#' @examples
#'
#' # Resolution range (up to 2 angstroms)
#' s <- seq(0,0.5,length.out=1000)
#'
#' # The chosen atom is nitrogen
#' atom_name <- "N"
#'
#' # Scattering factor
#' fN <- scafac(atom_name,s)
#'
#' # Plot
#'
#' @seealso
#' \code{\link{sfcoeffs}}
#'
#' @export
scafac <- function(atom_name,s) {
  # Pick tabulated coefficients
  lcoeffs <- sfcoeffs(atom_name)
  A <- lcoeffs$a
  B <- lcoeffs$b
  C <- lcoeffs$c

  # Change variables to adapt it to tabulated results
  s <- s/2 # Used in crystallography
  #s <- s/(4*pi) # Used in physical applications

  # Build function value as sum of gaussians
  ff <- A[1]*exp(-B[1]*s*s)+A[2]*exp(-B[2]*s*s)+
        A[3]*exp(-B[3]*s*s)+A[4]*exp(-B[4]*s*s)+C

  return(ff)
}


#' Calculate symmetry-equivalent positions
#'
#' Function to calculate all positions related by symmetry (including cell
#' centring) to the position of an atom in the asymmetric unit
#'
#' One atom has fractional coordinates \eqn{(x,y,z)}. The function first
#' finds the symmetry-equivalent \eqn{(x',y',z')} based on the operation,
#' \deqn{
#'   \left(\begin{array}{c}
#'     x' \\ y' \\ z'
#'   \end{array}\right) = R
#'   \left(\begin{array}{c}
#'     x' \\ y' \\ z'
#'   \end{array}\right) + T
#' }
#' where \eqn{R} is a 3 X 3 matrix representing the point group operation
#' and \eqn{T} a 3 X 1 vector representing the group translation. The
#' function then finds the copies \eqn{(x'',y'',z'')} of the symmetry-related
#' atoms \eqn{(x',y',z')} using the centring operator 3 X 1 vector \eqn{C}:
#' \deqn{
#'   \left(\begin{array}{c}
#'     x'' \\ y'' \\ z''
#'   \end{array}\right) =
#'   \left(\begin{array}{c}
#'     x' \\ y' \\ z'
#'   \end{array}\right) + C
#' }
#'
#' @param xyzf An n X 3 array or data frame whose rows are the fractional
#'             coordinates of the n atoms in the asymmetric unit.
#' @param aname A character string. The name of the atomic species.
#' @param B A vector of length n providing the thermal factor for the n
#'          atoms in the asymmetric unit.
#' @param Occ A vector of length n providing the occupancy numbers for the
#'            n atoms in the asymmetric unit.
#' @param SG A character string indicating the extended
#'           Hermann-Mauguin symbol for the space group.
#' @return set An integer equal to the specific setting
#'             corresponding to the given xHM symbol.
#' @param inau Logical variable. If TRUE, all atoms outside the unit cell
#'             are changed into atoms inside the unit cell by translational
#'             repetition (default is not to change expanded atoms).
#' @return An m X 6 data frame with column names x, y, z, atom, B, Occ,
#'         where m is the total number of atoms after the symmetry-expansion
#'         and the centring. The first three columns contain the three
#'         fractional coordinates, while the last three contain the atom name,
#'         the B factor and the occupancy.
#'
#' @examples
#'
#' # Create a "pretend" structure with five atoms in C 1 2 1
#' # Asymmetric unit is 0<=x<=1/2; 0<=y<1/2; 0<=z<1
#' xyzf <- matrix(nrow=5,ncol=3)
#' xyzf[1,] <- c(0.2,0.1,0.6)
#' xyzf[2,] <- c(0.15,0.15,0.55)
#' xyzf[3,] <- c(0.05,0.2,0.4)
#' xyzf[4,] <- c(0.25,0.15,0.56)
#' xyzf[5,] <- c(0.3,0.1,0.7)
#'
#' # Atom names
#' aname <- c("C","O","C","N","N")
#'
#' # Random B factors
#' B <- rnorm(5,mean=1,sd=0.2)
#'
#' # Occupancies
#' Occ <- rep(1,times=5)
#'
#' # Space group
#' SG <- "C 1 2 1"
#'
#' # Expansion
#' xyz <- expand_au(xyzf,aname,B,Occ,SG)
#' print(xyz) # The expanded structure should have 20 atoms
#'
#' @export
expand_au <- function(xyzf,aname,B,Occ,SG,inau=FALSE) {
  # Check input

  # Atoms in asymmetric unit
  dm <- dim(xyzf)
  if (length(dm) == 2) {
    if (dm[2] != 3) {
      msg <- "xyzf must be a n X 3 array of fractional atomic coordinates.\n"
      cat(msg)
      return(NULL)
    }
  } else {
    msg <- "xyzf must be a n X 3 array of fractional atomic coordinates.\n"
    cat(msg)
    return(NULL)
  }

  # Atom name
  if (!is.character(aname)) {
    msg <- "aname must be a characcter vector, the atom names.\n"
    cat(msg)
    return(NULL)
  }

  # B factors
  if (!is.numeric(B)) {
    msg <- "B must be a numeric vector of thermal factors.\n"
    cat(msg)
    return(NULL)
  }

  # Space group
  if (!is.character(SG)) {
    msg <- "SG must be a valid character symbol for an existing space group.\n"
    cat(msg)
    return(NULL)
  }
  if (is.null(findHM(SG))) {
    msg <- "SG must be a valid character symbol for an existing space group.\n"
    cat(msg)
    return(NULL)
  }

  # Occupancy
  if (!is.numeric(Occ)) {
    msg <- "Occ must be a numeric vector of atomic occupancies.\n"
    cat(msg)
    return(NULL)
  }

  # Check differences in length
  nxyzf <- length(xyzf[,2])
  naname <- length(aname)
  nB <- length(B)
  nOcc <- length(Occ)
  if (length(unique(c(nxyzf,naname,nB,nOcc))) != 1) {
    msg <- "Arrays xyzf, aname, B and Occ must have the same number of rows.\n"
    cat(msg)
    return(NULL)
  }

  # All checks passed

  # Operators related to symmetry
  lsymm <- syminfo_to_matrix_list(SG)

  # Number of symmetry operations
  nsym <- length(lsymm$PG)

  # Number of centring operations
  ncent <- length(lsymm$C)

  # Multiplicity
  nmult <- nsym*ncent

  # New array (data frame because of atom name)
  xyz <- matrix(rep(0,times=nxyzf*nsym*ncent),ncol=6,nrow=nxyzf*nmult)
  xyz <- as.data.frame(xyz)
  names(xyz) <- c("x","y","z","atom","B","Occ")

  # Fill with au
  xyz[1:nxyzf,1:3] <- xyzf
  xyz[1:nxyzf,4] <- aname
  xyz[1:nxyzf,5] <- B
  xyz[1:nxyzf,6] <- Occ

  # Expansion due to symmetry
  for (i in 2:nsym) {
    istart <- (i-1)*nxyzf+1
    iend <- i*nxyzf
    xyz[istart:iend,1:3] <- xyzf %*% t(lsymm$PG[[i]]) +
      matrix(rep(lsymm$T[[i]],times=nxyzf),ncol=3,byrow=TRUE)
    xyz[istart:iend,4] <- aname
    xyz[istart:iend,5] <- B
    xyz[istart:iend,6] <- Occ
  }
  ntot <- iend

  # Expansion due to centring (if available)
  ncent <- length(lsymm$C)
  if (ncent > 1) {
    for (i in 2:ncent) {
      istart <- (i-1)*ntot+1
      iend <- i*ntot
      tmp <- matrix(rep(lsymm$C[[i]],times=ntot),ncol=3,byrow=TRUE)
      xyz[istart:iend,1:3] <- xyz[1:ntot,1:3] + tmp[,1:3]
      xyz[istart:iend,4:6] <- xyz[1:ntot,4:6]
    }
  }

  # Atoms inside unit cell (if rquested)
  if (inau) {
    xyz[,1:3] <- xyz[,1:3] %% 1
  }

  return(xyz)
}


#' Calculation of structure factors
#'
#' Calculation of structure factors for any input reflection
#'
#' The theoretical formula in P 1 for the calculation of a structure factor
#' with Miller indices \eqn{(h,k,l)} is,
#' \deqn{
#'   F(h,k,l) = \sum_{j=1}^N Occ_j f_j(s)\exp(-B_j s^2/4)
#'   \exp(2\pi\mathrm{i}(hx_j+ky_j+lz_j))
#' }
#' where \eqn{B_j} is the B factor for atom j, \eqn{s=|\mathbf{h}|},
#' \eqn{f_j} is the scattering factor for atom j and \eqn{Occ_j} is the
#' occupancy for atom j.
#' When a symmetry is present, only the atoms in the asymmetric unit are
#' necessary as they are expanded into a P 1 structure using symmetry
#' operators. The specific symmetry will be reflected in amplitude and
#' phase of specific structure factors having precise values (e.g. 0, 90,
#' 180, etc) or being zero (systematic absences).
#'
#' @param hkl An m X 3 array or data frame of integers. Each row is a set
#'            of Miller indices, \eqn{(h,k,l)}.
#' @param xyzf An n X 3 array or data frame whose rows are the fractional
#'             coordinates of the n atoms in the asymmetric unit.
#' @param aname A character string. The name of the atomic species.
#' @param B A vector of length n providing the thermal factor for the n
#'          atoms in the asymmetric unit.
#' @param Occ A vector of length n providing the occupancy numbers for the
#'            n atoms in the asymmetric unit.
#' @param SG A character string indicating the extended
#'           Hermann-Mauguin symbol for the space group.
#' @param cpars A length-6 vector containing the unit cell parameters. The
#'              first three numbers are a, b, c, the last three are the
#'              angles \eqn{\alpha, \beta, \gamma} in degrees.
#' @return A named list with amplitudes (name "F") and phases in degrees
#'         (name "phi") of the structure factors calculated. Each component
#'         of the list is a vector of length equal to the number of input
#'         Miller indices.
#'
#' @examples
#'
#' # Random 5-atom structure in P 21 21 21
#' # Asymmetric unit is 0 < x < 1/4, 0 < y < 1, 0 < z < 1
#' xyzf <- matrix(c(runif(5,min=0,max=0.25),
#'                  runif(5,min=0,max=1),
#'                  runif(5,min=0,max=1)),ncol=3)
#'
#' # Atoms are C, N, O, F, Fe
#' aname <- c("C","N","O","F","Fe")
#'
#' # Random B factors around 1
#' B <- rnorm(5,mean=1,sd=0.5)
#'
#' # No atoms in special position
#' Occ <- rep(1,times=5)
#'
#' # Space group xHM symbol
#' SG <- "P 21 21 21"
#'
#' # Cell parameters
#' cpars <- c(20,30,40,90,90,90)
#'
#' # Create a few Miller indices
#' hkl <- expand.grid(h=-2:2,k=-2:2,l=-2:2)
#'
#' # Creation of structure factors
#' lSF <- strufac(hkl,xyzf,aname,B,Occ,SG,cpars)
#'
#' # Check systematic absences (give lSF$F = 0)
#' # For P 21 21 21 they are:
#' # (h,0,0),h=2n
#' # (0,k,0),k=2n
#' # (0,0,l),l=2n
#' fullidx <- 1:125
#' jdx <- sysabs(hkl,SG)
#' idx <- fullidx[-jdx]  # These are the systematic absences
#' print(hkl[idx,])
#' print(lSF$F[idx])
#'
#' # Now check reflections with (h,k,0). h even has phase 0,180
#' # and h odd has phase 90,270
#' idx <- which(hkl[,3] == 0)
#' print(lSF$phi[idx])
#'
#' @export
strufac <- function(hkl,xyzf,aname,B,Occ,SG,cpars) {
  # Check input

  # Miller indices
  dm <- dim(hkl)
  if (length(dm) == 2) {
    if (dm[2] != 3) {
      msg <- "hkl must be a n X 3 array of Miller indices.\n"
      cat(msg)
      return(NULL)
    }
  } else {
    msg <- "hkl must be a n X 3 array of Miller indices.\n"
    cat(msg)
    return(NULL)
  }
  nrefs <- length(hkl[,1])

  # Atoms in asymmetric unit
  if (!is.numeric(xyzf)) {
    msg <- "xyzf must be a n X 3 array of fractional atomic coordinates.\n"
    cat(msg)
    return(NULL)
  }
  dm <- dim(xyzf)
  if (length(dm) == 2) {
    if (dm[2] != 3) {
      msg <- "xyzf must be a n X 3 array of fractional atomic coordinates.\n"
      cat(msg)
      return(NULL)
    }
  } else {
    msg <- "xyzf must be a n X 3 array of fractional atomic coordinates.\n"
    cat(msg)
    return(NULL)
  }

  # B factors
  if (!is.numeric(B)) {
    msg <- "B must be a numeric vector of thermal factors.\n"
    cat(msg)
    return(NULL)
  }

  # Atom name
  if (!is.character(aname)) {
    msg <- "aname must be a characcter vector, the atom names.\n"
    cat(msg)
    return(NULL)
  }

  # Space group
  if (!is.character(SG)) {
    msg <- "SG must be a valid character symbol for an existing space group.\n"
    cat(msg)
    return(NULL)
  }
  if (is.null(findHM(SG))) {
    msg <- "SG must be a valid character symbol for an existing space group.\n"
    cat(msg)
    return(NULL)
  }

  # Occupancy
  if (!is.numeric(Occ)) {
    msg <- "Occ must be a numeric vector of atomic occupancies.\n"
    cat(msg)
    return(NULL)
  }

  # Check differences in length
  nxyzf <- length(xyzf[,2])
  naname <- length(aname)
  nB <- length(B)
  nOcc <- length(Occ)
  if (length(unique(c(nxyzf,naname,nB,nOcc))) != 1) {
    msg <- "Arrays xyzf, aname, B and Occ must have the same number of rows.\n"
    cat(msg)
    return(NULL)
  }

  # Cell parameters
  if (!is.numeric(cpars)) {
    msg <- "cpars must be a numeric vector with 6 cell parameters.\n"
    cat(msg)
    return(NULL)
  }
  if (length(cpars) != 6) {
    msg <- "cpars must be a numeric vector with 6 cell parameters.\n"
    cat(msg)
    return(NULL)
  }

  # All checks passed

  # Turn hkl into a matrix (for a later calculation)
  hkl <- as.matrix(hkl)

  # Expand asymmetric unit by symmetry
  xyz <- expand_au(xyzf,aname,B,Occ,SG)
  nmult <- length(xyz[,1])%/%length(aname)

  # Prepare quantities for SF (done my best to make it quick)
  parg <- 2*pi*hkl %*% t(xyz[,1:3])
  TRC <- cos(parg)
  TRS <- sin(parg)
  dd <- hkl_to_reso(hkl[,1],hkl[,2],hkl[,3],
                    cpars[1],cpars[2],cpars[3],
                    cpars[4],cpars[5],cpars[6])
  s <- 1/dd
  DW <- t(outer(xyz[,5],s,DebWal)) # This is relatively quick
  FJ <- t(outer(aname,s,scafacV))  # This is time consuming. Done just for au
  FJ <- matrix(rep(1,nmult),ncol=nmult) %x% FJ # Kronecker product is quick

  # Structure factors for an h,k,l are derived as row sums
  BigMc <- FJ * DW * TRC
  BigMs <- FJ * DW * TRS
  AH <- rowSums(BigMc)
  BH <- rowSums(BigMs)
  Fcplex <- complex(length.out=nrefs,real=AH,imaginary=BH)

  # Amplitudes and phases
  FF <- Mod(Fcplex)
  FF[FF < 1e-6] <- 0
  ff <- Arg(Fcplex)
  ff[FF < 1e-6] <- 0
  ff <- ff*180/pi
  ff <- ff %% 360

  return(list(F=FF,phi=ff))
}



#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#
## Auxiliary functions, only used here
#
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#
# Debye-Waller factor#
#
# exp(-B*s^2/4)
#
DebWal <- function(B,s) {
  rr <- exp(-B*s^2/4)

  return(rr)
}

#
# Vectorization of scafac
#
scafacV <- Vectorize(scafac)


