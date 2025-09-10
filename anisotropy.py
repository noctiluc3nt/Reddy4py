import numpy as np
import warnings

#' Invariant analysis of Reynolds stress tensor
#'
#'@description Invariant analysis of Reynolds stress tensor, calculation of Lumley and barycentric map coordinates and anisotropy
#'@param a11 R11 element of Reynolds stress tensor: \code{u_sd^2}   (scalar or vector)
#'@param a12 R12 element of Reynolds stress tensor: \code{cov(u,v)} (scalar or vector)
#'@param a13 R13 element of Reynolds stress tensor: \code{cov(u,w)} (scalar or vector)
#'@param a22 R22 element of Reynolds stress tensor: \code{v_sd^2}   (scalar or vector)
#'@param a23 R23 element of Reynolds stress tensor: \code{cov(v,w)} (scalar or vector)
#'@param a33 R33 element of Reynolds stress tensor: \code{w_sd^2}   (scalar or vector)
#'@param plot should the barycentric map be plotted? default \code{plot=FALSE}
#'
#'@return list containing \code{xb}, \code{yb}, \code{eta}, \code{xi}, all eigenvalues and eigenvectors (\code{eta}, \code{xi} are the coordinates of the Lumley triangle and \code{xb}, \code{yb} the coordinates of the barycentric map)
#'@export
#'
#'@examples
#'calc_anisotropy(1,0,0,1,0,1) #isotropic
#'calc_anisotropy(1,0,1,1,0,1) #anisotropic
#'
def calc_anisotropy(a11, a12, a13, a22, a23, a33, plot=False):
    #if len(set(map(len, [a11, a12, a13, a22, a23, a33]))) > 1:
    #    warnings.warn("The given elements of the Reynolds stress tensor are not of equal length.")
    n=1
    #unit matrix / Kronecker delta
    delta = np.eye(3)  
    #initialize
    eta = np.zeros(n)
    xi = np.zeros(n)
    xb = np.zeros(n)
    yb = np.zeros(n)
    evals = np.empty((n, 3))
    evecs = np.empty((n, 3, 3))
    #symmetry
    a21 = a12
    a31 = a13
    a32 = a23
    rey = np.array([[a11, a12, a13],
                    [a21, a22, a23],
                    [a31, a32, a33]])
    if not np.any(np.isnan(rey)):
        if np.sum(np.diag(rey)) != 0:
            B=rey/np.sum(np.diag(rey))-1/3*delta
            #diagonalize the anisotropy matrix (calculate the eigenvalues)
            evals_, evecs_ = np.linalg.eig(B) #invariant analysis
            Bdiag=np.diag(evals_)
            evs_sort=np.sort(evals_)[::-1]
            eta=(1/3*(evs_sort[0]**2+evs_sort[0]*evs_sort[1]+evs_sort[1]**2))**(1/2)
            dummyX=evs_sort[0]*evs_sort[1]*(evs_sort[0]+evs_sort[1])
            xi=-np.sign(dummyX)*(0.5*abs(dummyX))**(1/3)
            #barycentric map
            C1c = evs_sort[0]-evs_sort[1]
            C2c = 2*(evs_sort[1]-evs_sort[2])
            C3c = 3*(evs_sort[2])+1
            #coordinates
            xb = C1c+0.5*C3c
            yb = C3c*np.sqrt(3)/2
    out = {
        "xb": xb,
        "yb": yb,
        "eta": eta,
        "xi": xi,
        "eigenvalues": evals_,
        "eigenvectors": evecs_
    }
    if plot:
        plot_barycentric_map(out['xb'], out['yb'])
    return out
