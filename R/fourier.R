#
# This file is part of the cry package
#


#' Inverse Fourier synthesis
#'
#' From density to structure factors
#'
#' The electron density in the unit cell, \eqn{\rho(x,y,z)}, where
#' \eqn{(x,y,z)} are fractional coordinates, can be computed using
#' amplitudes and phases of the structure factors. The synthesis
#' formula is:
#'  \deqn{
#'   \rho(x,y,z) = \frac{1}{V}\sum_{h,k,l} F(h,k,l)
#'   \exp(-2\pi\mathrm{i}(hx+ky+lz))
#' }
#' where \eqn{F(h,k,l)} is the structure factor for the set of Miller indices
#' \eqn{h,k,l} and \eqn{V} is the unit cell volume. The inverse Fourier
#' synthesis yields the structure factors from the electron density.
#'
invfousynth3D <- function(rho,hkl) {
  # Number of grid points
  tmp <- dim(rho)
  nx <- tmp[1]
  ny <- tmp[2]
  nz <- tmp[3]
  N <- nx*ny*nz

  # Max frequencies
  H <- max(abs(hkl[,1]))
  K <- max(abs(hkl[,2]))
  L <- max(abs(hkl[,3]))

  # Nyquist limit
  if (H >= 0.5*nx | K >= 0.5*ny | L >= 0.5*nz) {
    stop("Grid too coarse. Decrease max Miller indices, or increase grid sampling.")
  }

  # All checks ok. Carry on.

  # Re-define max frequencies
  evorodd <- nx%%2
  if (evorodd == 0) {
    H <- as.integer(nx/2)
  } else {
    H <- as.integer((nx-1)/2)
  }
  evorodd <- ny%%2
  if (evorodd == 0) {
    K <- as.integer(ny/2)
  } else {
    K <- as.integer((ny-1)/2)
  }
  evorodd <- nz%%2
  if (evorodd == 0) {
    L <- as.integer(nz/2)
  } else {
    L <- as.integer((nz-1)/2)
  }

  # Create reference indices.
  # First (h,k,l) are transformed into (ii,jj,kk) with the
  # help of "sign", max h, max k, max l. Then a sequential
  # index is created using the standard convention of going
  # from cell [1,1,1] to cell [nx,ny,xz], column after column
  # and section after section.
  ii <- as.integer(hkl[,1] + 0.5*(1+sign(hkl[,1]+0.5))+
                   (1-sign(hkl[,1]+0.5))*(H+1))
  if (nx%%2 == 0) ii <- as.integer(hkl[,1] +
                        0.5*(1+sign(hkl[,1]+0.5))+
                        0.5*(1-sign(hkl[,1]+0.5))*(2*H+1))
  jj <- as.integer(hkl[,2] + 0.5*(1+sign(hkl[,2]+0.5))+
                     (1-sign(hkl[,2]+0.5))*(K+1))
  if (ny%%2 == 0) jj <- as.integer(hkl[,2] +
                        0.5*(1+sign(hkl[,2]+0.5))+
                        0.5*(1-sign(hkl[,2]+0.5))*(2*K+1))
  kk <- as.integer(hkl[,3] + 0.5*(1+sign(hkl[,3]+0.5))+
                     (1-sign(hkl[,3]+0.5))*(L+1))
  if (nz%%2 == 0) kk <- as.integer(hkl[,3] +
                        0.5*(1+sign(hkl[,3]+0.5))+
                        0.5*(1-sign(hkl[,3]+0.5))*(2*L+1))

  # Reference indices
  idx <- (kk-1)*nx*ny+(jj-1)*nx+ii
  #tmp <- c(hkl[,1],hkl[,2],hkl[,3])
  #print(matrix(c(tmp,ii,jj,kk,idx),ncol=7))

  # The actual FFT
  G <- fft(rho,inverse=TRUE)/N

  # Only interested in input indices
  FF <- G[idx]
  Fmod <- Mod(FF)
  Fphi <- Arg(FF)*180/pi
  jdx <- which(Fphi < 0)
  Fphi[jdx] <- Fphi[jdx]+360
  jdx <- which(Fphi >= 360)
  Fphi[jdx] <- Fphi[jdx]-360
  jdx <- which(abs(Im(FF)) < 0.000001 & abs(Re(FF)) >= 0.000001 & Re(FF) > 0)
  Fphi[jdx] <- 0.0
  jdx <- which(abs(Im(FF)) < 0.000001 & abs(Re(FF)) >= 0.000001 & Re(FF) < 0)
  Fphi[jdx] <- 180.0

  return(list(Fmod=Fmod,Fphi=Fphi,G=G,hkl=hkl))

  return(G)
}


fousynth3D <- function(hkl,Fmod,Fphi,nx,ny,nz,cpars) {
  # Expand in P 1
  #hklM <- rbind(hkl,-hkl)
  #FmodM <- matrix(c(Fmod,Fmod),ncol=1)  # Friedel's law
  #tmp <- -Fphi                          # Friedel's law
  #FphiM <- matrix(c(Fphi,tmp),ncol=1)
  #idx <- which(duplicated(hklM))
  #hkl <- hklM[-idx,]
  #Fmod <- FmodM[-idx]
  #Fphi <- FphiM[-idx]

  # Max frequencies
  H0 <- max(abs(hkl[,1]))
  K0 <- max(abs(hkl[,2]))
  L0 <- max(abs(hkl[,3]))

  # Nyquist limit
  if (H0 >= 0.5*nx | K0 >= 0.5*ny | L0 >= 0.5*nz) {
    stop("Grid too coarse. Decrease max Miller indices, or increase grid sampling.")
  }

  # Re-define max frequencies for the transform
  evorodd <- nx%%2
  if (evorodd == 0) {
    H <- as.integer(nx/2)
  } else {
    H <- as.integer((nx-1)/2)
  }
  evorodd <- ny%%2
  if (evorodd == 0) {
    K <- as.integer(ny/2)
  } else {
    K <- as.integer((ny-1)/2)
  }
  evorodd <- nz%%2
  if (evorodd == 0) {
    L <- as.integer(nz/2)
  } else {
    L <- as.integer((nz-1)/2)
  }

  # Create reference indices
  ii <- as.integer(hkl[,1] + 0.5*(1+sign(hkl[,1]+0.5))+
        (1-sign(hkl[,1]+0.5))*(H+1))
  if (nx%%2 == 0) ii <- as.integer(hkl[,1] +
                        0.5*(1+sign(hkl[,1]+0.5))+
                        0.5*(1-sign(hkl[,1]+0.5))*(2*H+1))
  jj <- as.integer(hkl[,2] + 0.5*(1+sign(hkl[,2]+0.5))+
                     (1-sign(hkl[,2]+0.5))*(K+1))
  if (ny%%2 == 0) jj <- as.integer(hkl[,2] +
                        0.5*(1+sign(hkl[,2]+0.5))+
                        0.5*(1-sign(hkl[,2]+0.5))*(2*K+1))
  kk <- as.integer(hkl[,3] + 0.5*(1+sign(hkl[,3]+0.5))+
                     (1-sign(hkl[,3]+0.5))*(L+1))
  if (nz%%2 == 0) kk <- as.integer(hkl[,3] +
                        0.5*(1+sign(hkl[,3]+0.5))+
                        0.5*(1-sign(hkl[,3]+0.5))*(2*L+1))

  # Reference indices
  idx <- (kk-1)*nx*ny+(jj-1)*nx+ii
  #tmp <- c(hkl[,1],hkl[,2],hkl[,3])
  #print(matrix(c(tmp,ii,jj,kk,idx),ncol=7))

  # Make structure factors as complex numbers
  FF <- complex(modulus=Fmod,argument=Fphi*pi/180,
                length.out=length(idx))

  # Volume of the unit cell
  ltmp <- lattice_stuff(cpars[1],cpars[2],cpars[3],
                        cpars[4],cpars[5],cpars[6])
  V <- ltmp[16]

  # Prepare full 3D array for fft
  G <- array(rep(0+0i,length=nx*ny*nz),c(nx,ny,nz))
  G[idx] <- FF/V  # Remember to divide by the unit cell volume

  # Density via fft
  rr <- Re(fft(G))

  # Return the appropriate grid
  x <- seq(0,1,length=nx+1)[1:nx]
  y <- seq(0,1,length=ny+1)[1:ny]
  z <- seq(0,1,length=nz+1)[1:nz]

  return(list(x=x,y=y,z=z,rho=rr,G=G))
}


#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#
## Auxiliary functions, only used here
#
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#
###
## Artificial 3D rho, for internal use only.
## Here just to test FFT. hkl must be a unique set
## of reflections
###
arti3Drho <- function(nx,ny,nz,hkl,Fmod,Fphi) {
  # 3D grid
  x <- seq(0,1,length=(nx+1))[1:nx]
  y <- seq(0,1,length=(ny+1))[1:ny]
  z <- seq(0,1,length=(nz+1))[1:nz]

  # Fill density matrix with a "trigonometric wave"
  N <- nx*ny*nz
  rho <- array(rep(Fmod[1],times=N),c(nx,ny,nz))
  for (i in 1:nx) {
    for (j in 1:ny) {
      for (k in 1:nz) {
        xx <- x[i]
        yy <- y[j]
        zz <- z[k]
        for (jdx in 2:length(hkl[,1])) {
          hh <- hkl[jdx,1]
          kk <- hkl[jdx,2]
          ll <- hkl[jdx,3]
          rho[i,j,k] <- rho[i,j,k]+2*Fmod[jdx]*
          cos(Fphi[jdx]*pi/180-2*pi*(hh*xx+kk*yy+ll*zz))
        }
      }
    }
  }

  return(list(x=x,y=y,z=z,rho=rho))
}
