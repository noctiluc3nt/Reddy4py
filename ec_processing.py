import numpy as np
import pandas as pd
import math

import os
os.chdir("/home/lauramack/clickhouse-db-data-processing/Reddy4py")

from auxillary import *
from constants import *
from diagnostics_meteorology import *
from diagnostics_turbulence import *


#os.chdir("./Reddy4py/")
#from auxillary import *
#from constants import *
#from diagnostics_meteorology import *
#from diagnostics_turbulence import *
#from ec_processing import *

#' Despiking
#'
#'@description Applies (up to) three despiking methods based on pre-defined thresholds
#'@param ts timeseries that shall be despiked
#'@param thresholds vector with two elements representing lower and upper bounds for despiking (pre-defined thresholds), \code{NA} means that the respective bound is not used
#'
#'@return despiked timeseries
#'@export
#'
#'@importFrom pracma detrend
#'
#'@examples
#'set.seed(5)
#'ts1=rnorm(100)
#'despiking(ts1,thresholds=c(-1,1))
#'
#'ts2=rexp(1000)
#'despiking(ts2)
#'
def despiking(ts,threshold_low,threshold_up):
    vec=np.where(ts<threshold_low,np.nan,ts)
    vec=np.where(ts>threshold_up,np.nan,ts)
    return ts

#' Count spikes
#'
#'@description Counts spikes in timeseries
#'
#'@param ts time series
#'@param thresholds vector with lower and upper threshold, e.g. c(0,10)
#'
#'@return number of spikes in timeseries (i.e. values lower than lower threshold and higher than upper threshold)
#'@export
#'
def count_spikes(ts,threshold_low,threshold_up):
    nr_spikes = np.sum(ts<threshold_low) + np.sum(ts>threshold_up)
    return nr_spikes

#' Amplitude resolution
#'
#'@description Gives amplitude resolution of time series (i.e. number of different values in time series)
#'
#'@param ts time series
#'
#'@return number of different values in time series
#'
def get_amplitude_resolution(ts):
    return len(set(ts))

#' Set everything smaller than machine epsilon to zero
#'
#'@description Calculates machine epsilon (machine-dependent) and sets everything smaller to exactly zero
#'@param vec vector/time series
#'
#'@return vector of same length, just all values smaller than machine epsilon are set to exactly zero
#'@export
#'
#'@examples
#'ts=c(1,0.1,1e-15,1e-16,1e-17,1e-18,1e-19)
#'ts=smaller_than_machine_epsilon(ts)
#'
def smaller_than_machine_epsilon(vec):
    epsilon=np.finfo(float).eps
    #2.22044604925e-16
    vec=np.where(vec<epsilon,0,vec)
    return vec

#' Double rotation
#'
#'@description Double rotation (i.e., sonic coordinate system will be aligned with streamlines)
#'@param u u-wind (levelled sonic)
#'@param v v-wind (levelled sonic)
#'@param w w-wind (levelled sonic)
#'
#'@return list containing the wind in a natural coordinate system (streamwise, crosswise, vertical) and the two rotation angles theta and phi
#'@export
#'
#'@examples
#'wind_rotated=rotate_double(4,3,1) #double rotation can be applied instantenously
#'
def rotate_double(u,v,w):
    #horizontal
    theta=math.atan2(np.nanmean(v),np.nanmean(u))
    u1=u*np.cos(theta) + v*np.sin(theta)
    v1=-u*np.sin(theta) + v*np.cos(theta)
    w1=w
    #vertical
    phi=math.atan2(np.mean(w1),np.mean(u1))
    u2=u1*np.cos(phi) + w1*np.sin(phi)
    v2=v1
    w2=-u1*np.sin(phi)+w1*np.cos(phi)
    return u2, smaller_than_machine_epsilon(v2), smaller_than_machine_epsilon(w2), (theta*180/np.pi+360)%360, (phi*180/np.pi+360)%360

#' Unit conversion of "parts-per" (molar mixing ratio) to density (for closed-path gas analyzer)
#'
#'@description Unit conversion of "parts-per" to density (for closed-path gas analyzer)
#'@param ppt measurement in parts per thousand [ppt]
#'@param T_mean temperature [K]
#'@param pres pressure [Pa]
#'@param e water vapor pressure [Pa]
#'@param gas which gas? can be either \code{H2O}, \code{CO2}, \code{CH4} (if \code{CO2}/\code{CH4} is selected, make sure that it's still in ppt and not ppm as usual)
#'
#'@return density of the gas [kg/m^3]
#'@export
#'
def ppt2rho(ppt,T_mean=288.15, pres = 101325, e = 0, gas="H2O"):
    Vd=Runiversal()*T_mean/(pres-e) #volume of dry air [m^3/mol]
    if gas == "H2O":
        return ppt/1000*M_H2O()/Vd
    elif gas == "CO2":
        return ppt/1000*M_CO2()/Vd
    elif gas == "CH4":
        return ppt/1000*M_CH4()/Vd

#' Conversion of molar concentration to density
#'
#'@description Conversion of molar concentration to density
#'@param c molar concentration in mol/m^3
#'@param gas which gas? can be either \code{H2O}, \code{CO2}, \code{CH4}
#'
#'@return density of the gas [kg/m^3]
#'@export
#'
def molarconcentration2density(c, gas="H2O"):
    if gas == "H2O":
        return(c*M_H2O())
    elif gas == "CO2":
        return(c*M_CO2())
    elif gas == "CH4":
        return(c*M_CH4())
    else:
        print("You selected a gas which is not available for the conversion here.")
     
#' Conversion of density to mixing ratio
#'
#'@description Conversion of density to mixing ratio
#'@param rho density [kg/m^3]
#'
#'@return mixing ratio of the gas [kg/kg]
#'@export
#'
def density2mixingratio(rho):
    return(rho/rhoAir())


#' Ts2T
#'
#'@description Converts sonic temperature Ts to temperature T
#'
#'@param Ts sonic temperature [K] (similar as virtual temperature)
#'@param q specific humidity [kg/kg]
#'
#'@return temperature [K]
#'@export
#'
def Ts2T(Ts,q):
    return Ts*(1+Rd()/Rv()*q)

#' SND and cross-wind correction of sensible heat flux
#'
#'@description SND and cross-wind correction of sensible heat flux: converts the buoyancy flux cov(w,Ts) (based on sonic temperature Ts) to sensible heat flux
#'@param Ts_mean sonic temperature [K] (averaged)
#'@param u_mean u-wind [m/s] (averaged)
#'@param v_mean v-wind [m/s] (averaged)
#'@param cov_uw cov(u,w) [m^2/s^2]
#'@param cov_vw cov(v,w) [m^2/s^2]
#'@param cov_wTs cov(Ts,w) [K*m/s] (buoyancy flux)
#'@param cov_qw cov(q,w) [kg/kg*m/s] (optional)
#'@param A constant used in cross-wind correction, default \code{A = 7/8} for CSAT3
#'@param B constant used in cross-wind correction, default \code{B = 7/8} for CSAT3
#'@param sos speed of sound [m/s], default \code{sos = csound()} corresponding to 343 m/s
#'
#'@return SND correction of sensible heat flux
#'@export
#'
def SNDcorrection(Ts_mean,u_mean,v_mean,cov_uw,cov_vw,cov_wTs,cov_qw=None,A=7/8,B=7/8):
    if cov_qw is not None:
        #cross-wind and SND correction
        return cov_wTs - 0.51*cov_qw + 2*Ts_mean/csound()**2*(A*u_mean*cov_uw + B*v_mean*cov_vw)
    else:
        #only cross-wind correction
        return cov_wTs + 2*Ts_mean/csound()**2*(A*u_mean*cov_uw + B*v_mean*cov_vw)


#' WPL correction for water vapor flux
#'
#'@description WPL correction: density correction for trace gas fluxes (i.e., converts volume- to mass-related quantity)
#'@param Ts_mean temperature [K] (sonic temperature or corrected temperature)
#'@param cov_wTs covariance cov(w,Ts) [m/s*K]
#'@param rhoq_mean measured water vapor density [kg/m^3]
#'@param cov_rhoqw covariance cov (w,rhow) [m/s*kg/m^3]
#'
#'@return WPL correction of respective flux
#'@export
#'
def WPLcorrectionH2O(cov_rhoqw,cov_wTs,Ts_mean,rhoq_mean):
    return((1+1.61*rhoq_mean/rhoAir())*(cov_rhoqw+rhoq_mean/Ts_mean*cov_wTs)) #with M_L/M_w = 1.61


#' WPL correction for CO2 flux
#'
#'@description WPL correction: density correction for trace gas fluxes (i.e., converts volume- to mass-related quantity)
#'@param Ts_mean temperature [K] (sonic temperature or corrected temperature)
#'@param cov_wTs covariance cov(w,Ts) [m/s*K]
#'@param rhoq_mean measured water vapor density [kg/m^3]
#'@param cov_rhoqw covariance cov (w,rhow) [m/s*kg/m^3]
#'
#'@return WPL correction of respective flux
#'@export
#'
def WPLcorrectionCO2(cov_rhocw,cov_rhoqw,cov_wTs,Ts_mean,rhoq_mean,rhoc_mean):
    return(cov_rhocw+1.61*rhoc_mean/rhoAir()*cov_rhoqw+(1+1.61*rhoq_mean/rhoAir())*rhoc_mean/Ts_mean*cov_wTs) #with M_L/M_w = 1.61



### speed of sound to sonic temperature ###
#' Converts speed of sound (sos) to sonic temperature
#'
#'@description Converts speed of sound (sos) to sonic temperature
#'@param sos speed of sound [m/s]
#'
#'@return sonic temperature (virtual temperature) [K]
#'@export
#' 
def sos2Ts(sos):
	return(sos**2/(cpcv()*Rd()))


#' Vertical Velocity Flag
#'
#'@description Vertical velocity flag according to Mauder et al., 2013: After rotation the vertical velocity should vanish, this flag flags high remaining vertical velocities.
#'@param w vertical velocity
#'@param thresholds_w vector containing 2 elements to distinguish between flag=0 and flag=1, as well as flag=1 and flag=2, default: \code{c(0.1,0.15)}
#'
#'@return vertical velocity flags (0: in full agreement with the criterion ... 2: does not fulfill the criterion)
#'@export
#'
#'@examples
#'flag_w(0.01)
#'
def flag_w(w,thresholds_w=[0.1,0.15]):
    if (len(thresholds_w)!=2):
        print("thresholds_w has to be a vector of length 2.")
    w=abs(w)
    flag=np.where(w<thresholds_w[0],0,np.where(w<thresholds_w[1],1,2))
    return(flag)



#' Integral Turbulence Characteristics Flag 
#'
#'@description Integral Turbulence Characteristics Flag: Tests the consistency with Monin-Obukhov similarity theory using the scaling functions from Panofsky and Dutton, 1984.
#'@param w_sd standard deviation of vertical velocity
#'@param ustar friction velocity
#'@param zeta stability parameter \code{zeta = z/L}
#'@param thresholds_most vector containing 2 elements to distinguish between flag=0 and flag=1, as well as flag=1 and flag=2, default: \code{c(0.3,0.8)}
#'
#'@return integral turbulence characteristics flags (0: in full agreement with the criterion ... 2: does not fulfill the criterion)
#'@export
#'
#'@examples
#'itc_flag=flag_most(0.2,0.4,-0.3)
#'
def flag_most(w_sd,ustar,zeta,thresholds_most=[0.3,0.8]):
    if (len(thresholds_most)!=2):
        print("thresholds_most has to be a vector of length 2.")
    parameterized=1.3*(1-2*abs(zeta))**(1/3) #w_sd/ustar parametrized according to scaling function based on zeta
    itc=abs((w_sd/ustar-parameterized)/parameterized)
    flag=np.where(itc<thresholds_most[0],0,np.where(itc<thresholds_most[1],1,2))
    return(flag)



### fluxes (unit conversion) ###
#' Converts cov(w,T) to sensible heat flux SH
#'
#'@description Converts cov(T,w) to sensible heat flux SH
#'@param cov_wT covariance cov(w,T) [K m/s]
#'@param rho density of air [kg/m^3] (optional)
#'
#'@return sensible heat flux [W/m^2]
#'@export
#'
def cov2sh(cov_wT,rho=None):
	if rho is None:
		rho = rhoAir()
	return(rho*cp()*cov_wT)


#' Converts cov(w,q) to latent heat flux LH
#'
#'@description Converts cov(w,q) to latent heat flux LH
#'@param cov_wq covariance cov(w,q) [m/s]
#'@param rho density of air [kg/m^3] (optional)
#'
#'@return latent heat flux [W/m^2]
#'@export
#'
def cov2lh(cov_wq,rho=None):
	if (rho is None):
		rho = rhoAir()
	return(rho*Lv()*cov_wq)


#' Converts cov(co2,w) to CO2 flux
#'
#'@description Converts cov(co2,w) to CO2 flux
#'@param cov_co2w covariance cov(co2,w) [m/s]
#'@param rho density of air [kg/m^3] (optional)
#'
#'@return CO2 flux [kg/(m^2*s)]
#'@export
#'
def cov2cf(cov_co2w,rho=None):
    if rho is None:
        rho = rhoAir()
    return(rho*cov_co2w)
