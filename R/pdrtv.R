#
# This file is part of the cry package
#
# Functions connected to reflections data.

#' Reads and output a CIF file for powder diffraction
#'
#' @param filename A character string. The path to a valid CIF file.
#' @param message A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the cif file.
#' @return A named list. Each name correspond to a valid field in the powder
#'    diffraction Rietveld processed CIF.
#'
#' @examples
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"e-65-00i60-Isup2.rtv")
#' lCIF <- readpd_rtv(filename)
#' print(names(lCIF))
#' print(lCIF$INTRO$CELL)
#' print(lCIF$INTRO$HALL)
#' print(lCIF$INTRO$HM)
#' print(lCIF$REFL)
#'
#' @export
readpd_rtv <- function(filename,message=FALSE){
  f <- file(filename)
  lcif <- readLines(f,warn=FALSE)
  #c_lcif <- lcif[!grepl('^#|^.*#', lcif)]
  l_list <- grep("loop_",lcif)
  l_list1 <- l_find(l_list, length(lcif))
  mat<-zoo::rollapply(l_list1, 2,stack)
  ch <- apply(mat, 1, function(x) lcif[(x[1]+1):(x[2]-1)])
  diffraction <- lapply(ch, ucheck, pattern="_pd_proc_point_id")
  reflection <- lapply(ch, ucheck, pattern="_refln_index_h")

  intro <- r_summ(lcif)
  diffractions <- if (is.na(nanona(diffraction)) == FALSE) clean(r_diffractions(nanona(diffraction))) else NULL
  reflections <- if (is.na(nanona(reflection)) == FALSE) clean(r_peakreflns(nanona(reflection))) else NULL
  CIF = list(HEADER=intro,DIFF=diffractions,REFL_PEAK=reflections)
  close(f)
  nrefs <- length(reflections$VAL$F_meas_au)
  fmeas <- as.numeric(reflections$VAL$F_meas_au)
  msg <- c("\n")
  if (message) {
    msg <- c(msg,sprintf("File %s read successfully.\n",filename))
    msg2 <- sprintf("There are %d reflections in this file.\n",nrefs)
    msg <- c(msg,msg2)
    msg <- c(msg,"Here is a summary of the observations:\n")
    msg <- c(msg,"\n")
    cat(msg)
    print(summary(fmeas))
  }

  return(CIF)
}



### accessory functions ####

l_find <- function(a,n){
  if (length(a) > 1){
    a1 <- append(a,(n-1))
  } else
  { a1 <- c(2,a,(n-1))
  return(a1)
  }
}

stack<-function(x){
  j <- c(x[1],x[2])
  return(j)
}

ucheck <- function(x,pattern){
  r <- unlist(x)
  if (length(grep(pattern,r))>0){
    piece <- r
  } else
  { piece <- NA
  return(piece)
  }
}

ansnull <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NULL
  } else
  { out <- x[!is.na(x)]
  return(out)
  }
}

nanona <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NA
  } else
  { out <- x[!is.na(x)]
  return(out)
  }
}

recheck <- function(r1){
  r2 <- gsub("[:):]","",gsub("[:(:]",",",r1))
  r3 <- as.numeric((strsplit(r2, ",")[[1]])[1])
  return(r3)
}

check <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "\\s+")[[1]])[2] else NA
  r2 <- if (length(grep("[:(:]",r1,value = TRUE)>0) == TRUE) recheck(r1) else r1
  return(r2)
}

check1 <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "'")[[1]])[2] else NA
  return(r1)
}

clean1 <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NULL
  } else
  { out <- nc_type(as.data.frame(x))
  return(out)
  }
}


clean <- function(x){
  co1 <- data.frame(gsub ("[()]","",as.matrix(x),perl=T),stringsAsFactors = FALSE)
  ref <- data.frame(gsub("(?<!\\))(?:\\w+|[^()])(?!\\))","",as.matrix(x),perl=T))
  ref1 <- data.frame(gsub("[()]","",as.matrix(ref),perl=T),stringsAsFactors = FALSE)
  ref1[ref1==""]<-NA
  ref2 <- clean1(ref1)
  col1 <- nc_type(co1)
  return(list(VAL=col1,STD=ref2))
}

reap1 <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NULL
  } else
  { out <- as.numeric(x)
  return(out)
  }
}

reap <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "\\s+")[[1]])[2] else NA
  v <- as.numeric(gsub ("[()]","",as.matrix(r1),perl=T))
  s <- gsub("(?<!\\))(?:\\w+|[^()])(?!\\))","",as.matrix(r1),perl=T)
  s1 <- gsub("[()]","",as.matrix(s),perl=T)
  s2 <- reap1(s1)
  return(list(VAL=v,STD=s2))
}

nc_type <- function(data){
  count <- as.numeric(ncol(data))
  if (isTRUE(count > 2)) {
    data[] <- lapply(data, function(x) numas(x))
    out <- data
  } else if (count == 2){
    l1_data <- list(data$VAL)
    l_1 <- lapply(l1_data[[1]], function(x) numas(x))
    l2_data <- list(data$KEY)
    l2 <- c(gsub("\\[|\\]" ,"",l2_data[[1]]))
    names(l_1) <- c(l2)
    out <- l_1
    return(out)
  }
}

numas <- function(x){
  data <- x
  out <- (suppressWarnings(as.numeric(data)))
  if (all(is.na(out))== FALSE) {
    out1 <- out
  } else {
    out1 <- as.character(data)
  }
  return(out1)
}

r_diffractions <- function (x){
  data <- unlist(x)
  data <- data[!grepl("_number_of_points", data)]
  data <- data[data !=" "]
  l_l <- c(grep("_pd",data,value=TRUE))
  m <- length(l_l)
  n <- length(data)
  data1 <- data[m+1:n]
  data1 <- nanona(data1)
  data2 <- scan(text=data1, what='character', quiet=TRUE)
  o <- length(data2)
  data3 <- list(data2)
  list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
  res <- do.call(rbind, list_all)
  res <- as.data.frame(res)
  l_l <- c(gsub(" ","",l_l))
  colnames(res) <- c(gsub("_pd_","",l_l))
  return(res)
}


r_peakreflns <- function (x){
  data <- unlist(x)
  data1 <- data[!grepl('^#|^.*#', data)]
  data1 <- data1[data1 !=" "]
  l_l <- c(grep("_refln_",data1,value=TRUE))
  m <- length(l_l)
  n <- length(data1)
  data2 <- data1[m+1:n]
  data2 <- nanona(data2)
  data3 <- scan(text=data2, what='character', quiet=TRUE)
  o <- length(data3)
  data4 <- list(data3)
  list_all <- split(data4[[1]], rep(1:(length(data3)/length(l_l)), each = m))
  res <- do.call(rbind, list_all)
  res <- as.data.frame(res)
  l_l <- c(gsub(" ","",l_l))
  colnames(res) <- c(gsub("_refln_","",l_l))
  return(res)
}

r_summ <- function(x){
  data <- unlist(x)
  id <- check("_pd_block_id",data)
  pdm <- check("_pd_calc_method",data)
  mtmi <- reap("_meas_2theta_range_min",data)
  mtma <- reap("_meas_2theta_range_max",data)
  mti <- reap("_meas_2theta_range_inc",data)
  ptmi <- reap("_proc_2theta_range_min",data)
  ptma <- reap("_proc_2theta_range_max",data)
  pti <- reap("_proc_2theta_range_inc",data)
  rp <- check("_diffrn_radiation_probe ",data)
  wl <-reap("_diffrn_radiation_wavelength ",data)
  rf <- reap("_proc_ls_prof_R_factor",data)
  wr <- reap("_proc_ls_prof_wR_factor",data)
  e_wr <- reap("_proc_ls_prof_wR_expected",data)
  Theta_r <- list(Meas_min=mtmi,Meas_max=mtma,Meas_inc=mti,Proc_min=ptmi,Proc_max=ptma,Proc_inc=pti)
  R_fit <- list(Rp=rf,Rwp=wr,Rexp=e_wr)
  summ_c <- list(ENTRY=id,THETA_RANGE=Theta_r,PROBE_TYPE=rp,WAVELENGTH=wl,R_PRO_FIT=R_fit)
  return(summ_c)
}
