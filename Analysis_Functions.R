# This file contains custom functions for peak detection and calculation of electrophysiological params


# Custom peak detection function

peakdet <- function(v, delta, x = NULL)
{
  maxtab <- NULL
  mintab <- NULL
  
  if (is.null(x))
  {
    x <- seq_along(v)
  }
  
  if (length(v) != length(x))
  {
    stop("Input vectors v and x must have the same length")
  }
  
  if (!is.numeric(delta))
  {
    stop("Input argument delta must be numeric")
  }
  
  if (delta <= 0)
  {
    stop("Input argument delta must be positive")
  }
  
  mn <- Inf
  mx <- -Inf
  
  mnpos <- NA
  mxpos <- NA
  
  lookformax <- TRUE
  
  for(i in seq_along(v))
  {
    this <- v[i]
    
    if (this > mx)
    {
      mx <- this
      mxpos <- x[i]
    }
    
    if (this < mn)
    {
      mn <- this
      mnpos <- x[i]
    }
    
    if (lookformax)
    {
      if (this < mx - delta)
      {
        maxtab <- rbind(maxtab, data.frame(pos = mxpos, val = mx))
        
        mn <- this
        mnpos <- x[i]
        
        lookformax <- FALSE
      }
    }
    else
    {
      if (this > mn + delta)
      {
        mintab <- rbind(mintab, data.frame(pos = mnpos, val = mn))
        
        mx <- this
        mxpos <- x[i]
        
        lookformax <- TRUE
      }
    }
  }
  
  list(maxtab = maxtab, mintab = mintab)
}