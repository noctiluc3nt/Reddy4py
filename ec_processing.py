def despiking(vec,threshold_low,threshold_up):
    vec=np.where(vec<threshold_low,np.nan,vec)
    vec=np.where(vec>threshold_up,np.nan,vec)
    return vec

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
    return u2, smaller_than_machine_epsilon(v2), smaller_than_machine_epsilon(w2)

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

