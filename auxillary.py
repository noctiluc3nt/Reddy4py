import numpy as np

#' Calculates circular mean
#'
#'@description calculates circular mean
#'@param x input vector, e.g. wind directions [degree]
#'@param na.rm should NA values be removed? default \code{TRUE}
#'@return circular mean of x values
#'@export
#'
#'@examples
#'wd=c(280,90)
#'calc_circular_mean(wd)
#'
def calc_circular_mean(x,nan_rm=True):
    x=x*np.pi/180
    if nan_rm is True:
        return((math.atan2(np.nansum(np.sin(x)),np.nansum(np.cos(x)))*180/np.pi)%360)
    else:
        return((math.atan2(np.sum(np.sin(x)),sum(np.cos(x)))*180/np.pi)%360)

