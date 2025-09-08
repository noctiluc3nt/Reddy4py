import numpy as np

def despiking(ts,threshold_low,threshold_up):
    vec=np.where(ts<threshold_low,np.nan,ts)
    vec=np.where(ts>threshold_up,np.nan,ts)
    return ts

def count_spikes(ts,threshold_low,threshold_up):
    nr_spikes = np.sum(ts<threshold_low) + np.sum(ts>threshold_up)
    return nr_spikes

def get_amplitude_resolution(ts):
    return len(set(ts))

def smaller_than_machine_epsilon(vec):
    epsilon=np.finfo(float).eps
    #2.22044604925e-16
    vec=np.where(vec<epsilon,0,vec)
    return vec

def rotate_double(u,v,w):
    #horizontal
    theta=math.atan2(np.mean(v),np.mean(u))
    u1=u*np.cos(theta) + v*np.sin(theta)
    v1=-u*np.sin(theta) + v*np.cos(theta)
    w1=w
    #vertical
    phi=math.atan2(np.mean(w1),np.mean(u1))
    u2=u1*np.cos(phi) + w1*np.sin(phi)
    v2=v1
    w2=-u1*np.sin(phi)+w1*np.cos(phi)
    return u2, smaller_than_machine_epsilon(v2), smaller_than_machine_epsilon(w2), (theta*180/np.pi+360)%360, (phi*180/np.pi+360)%360

def ppt2rho(ppt,T_mean=288.15, pres = 101325, e = 0, gas="H2O"):
    Vd=Runiversal()*T_mean/(pres-e) #volume of dry air [m^3/mol]
    if gas == "H2O":
        return ppt/1000*M_H2O()/Vd
    elif gas == "CO2":
        return ppt/1000*M_CO2()/Vd
    elif gas == "CH4":
        return ppt/1000*M_CH4()/Vd

def Ts2T(Ts,q):
    return Ts*(1+Rd()/Rv())*q

def SNDcorrection(Ts_mean,u_mean,v_mean,cov_uw,cov_vw,cov_wTs,cov_qw,A=7/8,B=7/8):
    if cov_qw is not None:
        #cross-wind and SND correction
        cov_wTs - 0.51*cov_qw + 2*Ts_mean/csound()**2*(A*u_mean*cov_uw + B*v_mean*cov_vw)
    else:
        #only cross-wind correction
        return cov_wTs + 2*Ts_mean/csound()**2*(A*u_mean*cov_uw + B*v_mean*cov_vw)
    
def WPLcorrection(Ts_mean,q_mean,cov_wTs,rhow_mean,cov_wrhow,rhoc_mean=None,cov_wrhoc=None):
    if rho_c is None: #water vapor flux
        return((1+1.61*q_mean)*(cov_wrhow+rhow_mean/Ts_mean*cov_wTs)) #with M_L/M_w = 1.61
    else: #other trace gas flux
        return(cov_wrhoc+1.61*rhoc_mean/rhow_mean*cov_wrhow+(1+1.61*q_mean)*rhoc_mean/Ts_mean*cov_wTs)
    
def sos2Ts(sos):
	return(sos^2/(cpcv()*Rd()))


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
