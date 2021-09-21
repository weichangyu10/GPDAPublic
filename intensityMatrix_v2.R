intensityMatrix_v2 <- function(spectra,t) {
#intensityMatrix <- function(peaks, spectra) {
  
#  ## test arguments
#  .stopIfNotIsMassPeaksList(peaks)
  
  #m <- .as.matrix.MassObjectList(spectra)
  N = length(spectra)
  TP = length(t)
  m = matrix(0,N, TP)
  colnames(m) = t
  
  for (i in 1:N) {
       approxSpectra <- approxfun(x = mass(spectra[[i]]), y = intensity(spectra[[i]]), yleft=0L, yright=0L)
       m[i,] <- approxSpectra(t)
     }
  
  # ## lookup corresponding intensity values in spectra for missing peaks
  # if (!missing(spectra)) {
  #   .stopIfNotIsMassSpectrumList(spectra)
  #   
  #   if (length(peaks) != length(spectra)) {
  #     stop("Incompatible number of spectra!")
  #   }
  #   
  #   isNa <- is.na(m)
  #   uniqueMass <- as.double(colnames(m))
  #   
  #   approxSpectra <- lapply(spectra, approxfun, yleft=0L, yright=0L)
  #   
  #   for (i in seq_along(approxSpectra)) {
  #     m[i, isNa[i, ]] <- approxSpectra[[i]](uniqueMass[isNa[i, ]])
  #   }
  # }
  
  m
}