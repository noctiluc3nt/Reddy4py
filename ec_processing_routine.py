import numpy as np
import pandas as pd
import os
from scipy import signal
import logging

os.chdir("/home/lauramack/flux-processing/")
from Reddy4py import *

#from Reddy4py.anisotropy import *
#from Reddy4py.auxillary import *
#from Reddy4py.constants import *
#from Reddy4py.diagnostics_meteorology import *
#from Reddy4py.diagnostics_turbulence import *
#from Reddy4py.ec_processing import *


def ec_processing_rt(dat,TIMESTAMP,convert_sos2Ts=True,
                     measurement_height=4.4,do_log=True):
    """
    eddy-covariance post-processing of high-frequency data from Finse 

    @input numpy array with columns: u-v-w-Ts-h2o-co2(-ch4) (high-frequency data from table {station}_HFData)
    @return post-processed, quality-controlled, averaged fluxes and co (for inserting back to table ch_{station}postproc)
    
    h2o, co2, ch4: uses molar concentrations [mol/m^3] and converted them to density [kg/m^3] by multiplying with molar mass
    
    """
    dat = np.array(dat,dtype=float,copy=True)
    if dat.shape[1] not in (6,7):
        raise ValueError("expected 6 columns (u, v, w, SOS/Ts, H2O, CO2) or 7 columns with CH4")
    #units
    if convert_sos2Ts:
        dat[:,3]=sos2Ts(dat[:,3])
    #dat[:,4]=ppt2rho(dat[:,4],gas="H2O") #ppt
    #dat[:,5]=ppt2rho(dat[:,5]/1000,gas="CO2") #ppm originally
    dat[:,4]=molarconcentration2density(dat[:,4]/1000,gas="H2O") #mmol/m^3 -> kg/m^3
    dat[:,5]=molarconcentration2density(dat[:,5]/1000,gas="CO2")
    #despiking + spike counting
    nr_spikes_u=count_spikes(dat[:,0],-15,15)
    dat[:,0]=despiking(dat[:,0],-15,15)
    nr_spikes_v=count_spikes(dat[:,1],-15,15)
    dat[:,1]=despiking(dat[:,1],-15,15)
    nr_spikes_w=count_spikes(dat[:,2],-5,5)
    dat[:,2]=despiking(dat[:,2],-5,5)
    nr_spikes_Ts=count_spikes(dat[:,3],230,300)
    dat[:,3]=despiking(dat[:,3],230,300)
    nr_spikes_h2o=count_spikes(dat[:,4],0,0.05)
    dat[:,4]=despiking(dat[:,4],0,0.05)
    nr_spikes_co2=count_spikes(dat[:,5],0,0.005)
    dat[:,5]=despiking(dat[:,5],0,0.005)
    #amplitude resolution
    ampl_res_u=get_amplitude_resolution(dat[:,0])
    ampl_res_v=get_amplitude_resolution(dat[:,1])
    ampl_res_w=get_amplitude_resolution(dat[:,2])
    ampl_res_Ts=get_amplitude_resolution(dat[:,3])
    ampl_res_h2o=get_amplitude_resolution(dat[:,4])
    ampl_res_co2=get_amplitude_resolution(dat[:,5])
    #wind direction + mean
    wd_mean=calc_circular_mean(calc_windDirection(dat[:,0],dat[:,1]))
    ws_mean=np.nanmean(calc_windspeed(dat[:,0],dat[:,1]))
    #double rotation
    wind_rotated=rotate_double(dat[:,0],dat[:,1],dat[:,2])
    dat[:,0],dat[:,1],dat[:,2]=wind_rotated[0:3]
    dr_rot1=wind_rotated[3]
    dr_rot2=wind_rotated[4]
    #if methane
    do_methane=(dat.shape[1]>6)
    if do_methane: 
        dat[:,6]=molarconcentration2density(dat[:,6]/1000,gas="CH4")
        nr_spikes_ch4=count_spikes(dat[:,6],0,3)
        dat[:,6]=despiking(dat[:,6],0,3)
        ampl_res_ch4=get_amplitude_resolution(dat[:,6])
        mask = ~np.isnan(dat[:,2]) & ~np.isnan(dat[:,6])
        #cov_ch4w=np.cov(dat[mask,2],dat[mask,6],bias=True)[0,1]
        cov_ch4w=calc_cov(dat[:,2],dat[:,6],maxlag=120)["cov_maximized"]
        #eddy stats use pandas quickly
        dat=pd.DataFrame(dat)
        #dat=dat.asfreq(freq='200L')
        u_mean,v_mean,w_mean,Ts_mean,h2o_mean,co2_mean,ch4_mean=dat.apply(np.nanmean)
        u_sd,v_sd,w_sd,Ts_sd,h2o_sd,co2_sd,ch4_sd=dat.apply(np.nanstd)
    else:
        #eddy stats use pandas quickly
        dat=pd.DataFrame(dat)
        #dat=dat.asfreq(freq='200L')
        u_mean,v_mean,w_mean,Ts_mean,h2o_mean,co2_mean=dat.apply(np.nanmean)
        u_sd,v_sd,w_sd,Ts_sd,h2o_sd,co2_sd=dat.apply(np.nanstd)
    #back to numpy
    dat=np.array(dat)
    q_mean=np.nanmean(dat[:,4])/1.225 #kg/kg
    T_mean=Ts2T(Ts_mean,q_mean)
    #mask = ~np.isnan(dat[:,2]) & ~np.isnan(dat[:,0])
    #cov_uw=np.cov(dat[mask,0],dat[mask,2],bias=True)[0,1]
    #mask = ~np.isnan(dat[:,0]) & ~np.isnan(dat[:,1])
    #cov_uv=np.cov(dat[mask,0],dat[mask,1],bias=True)[0,1]
    #mask = ~np.isnan(dat[:,2]) & ~np.isnan(dat[:,1])
    #cov_vw=np.cov(dat[mask,1],dat[mask,2],bias=True)[0,1]
    #mask = ~np.isnan(dat[:,2]) & ~np.isnan(dat[:,3])
    #cov_wTs=np.cov(dat[mask,2],dat[mask,3],bias=True)[0,1]
    #mask = ~np.isnan(dat[:,2]) & ~np.isnan(dat[:,4])
    #cov_h2ow=np.cov(dat[mask,2],dat[mask,4],bias=True)[0,1]
    #mask = ~np.isnan(dat[:,2]) & ~np.isnan(dat[:,5])
    #cov_co2w=np.cov(dat[mask,2],dat[mask,5],bias=True)[0,1]
    #with maximize cross-correlation
    cov_uw=calc_cov_pairwise_complete(dat[:,0],dat[:,2])
    cov_uv=calc_cov_pairwise_complete(dat[:,0],dat[:,1])
    cov_vw=calc_cov_pairwise_complete(dat[:,1],dat[:,2])
    cov_wTs=calc_cov_pairwise_complete(dat[:,2],dat[:,3])
    cov_h2ow=calc_cov(dat[:,2],dat[:,4],maxlag=120)["cov_maximized"]
    cov_co2w=calc_cov(dat[:,2],dat[:,5],maxlag=120)["cov_maximized"]
    #snd correction
    #cov_wT=SNDcorrection(Ts_mean,u_mean,v_mean,cov_uw,cov_vw,cov_wTs)
    cov_wT=cov_wTs
    #wpl correction
    #cov_h2ow=WPLcorrectionH2O(cov_h2ow,cov_wTs,Ts_mean,h2o_mean)
    #cov_co2w=WPLcorrectionCO2(cov_co2w,cov_h2ow,cov_wTs,Ts_mean,h2o_mean,co2_mean)
    #some more internal turbulence quantities
    sh=cov2sh(cov_wT)
    lh=cov2lh(cov_h2ow)
    et=lh2et(lh)
    br=calc_br(sh,lh)
    cf=cov2cf(cov_co2w)*1000 #g/(m^2*s)
    tke=calc_tke(u_sd,v_sd,w_sd)
    ustar=calc_ustar(cov_uw,cov_vw)
    z0=ustar2z0(ustar)
    dshear=calc_dshear(cov_uw,cov_vw)
    obukhov_length=calc_L(ustar,T_mean,cov_wT)
    zeta=calc_zeta(measurement_height,obukhov_length)
    flux_intermittency=calc_flux_intermittency(dat[:,3]) #intermittency in temperature
    #quality flags
    qf_w=flag_w(w_mean)
    qf_most=flag_most(w_sd,ustar,zeta)
    qf_stationarity=2
    try:
        qf_stationarity=flag_stationarity(dat[:,2],dat[:,3]) #stationarity flag for w'Ts'
    except:
        print("An error in stationary flag calculation occurred.")
        if do_log:
            logging.exception("Stationarity flag calculation failed for %s", TIMESTAMP)
    qf_all=max(qf_w,qf_most,qf_stationarity)
    #anisotropy
    xb=0
    yb=0
    try:
        aniso=calc_anisotropy(u_sd**2,cov_uv,cov_uw,v_sd**2,cov_vw,w_sd**2)
        xb=aniso['xb']
        yb=aniso['yb']
    except:
        print("An error in anisotropy calculation occurred.")
        if do_log:
            logging.exception("Anisotropy calculation failed for %s", TIMESTAMP)
        xb=np.nan
        yb=np.nan
    ### return
    dr_rot1 = np.int32(-1) if np.isnan(dr_rot1) else np.int32(dr_rot1)
    dr_rot2 = np.int32(-1) if np.isnan(dr_rot2) else np.int32(dr_rot2)
    row=pd.DataFrame([TIMESTAMP,u_mean,v_mean,w_mean,Ts_mean,T_mean,h2o_mean,co2_mean,
                     u_sd,v_sd,w_sd,Ts_sd,h2o_sd,co2_sd,
                     cov_uw,cov_vw,cov_uv,wd_mean,ws_mean,
                     ustar,tke,dshear,z0,
                     sh,lh,et,br,cf,
                     obukhov_length,zeta,xb,yb,flux_intermittency,
                     dr_rot1,dr_rot2,nr_spikes_u,nr_spikes_v,nr_spikes_w,
                     nr_spikes_Ts,nr_spikes_h2o,nr_spikes_co2,
                     ampl_res_u,ampl_res_v,ampl_res_w,ampl_res_Ts,ampl_res_h2o,ampl_res_co2,
                     qf_most,qf_stationarity,qf_w,qf_all]).T
    row.columns=["time","u_mean","v_mean","w_mean","Ts_mean","T_mean","h2o_mean","co2_mean",
                     "u_sd","v_sd","w_sd","Ts_sd","h2o_sd","co2_sd",
                     "cov_uw","cov_vw","cov_uv","wd_mean","ws_mean",
                     "ustar","tke","dshear","z0",
                     "SH","LH","ET","BR","CF",
                     "L","zeta","xb","yb","flux_intermittency",
                     "dr_rot1","dr_rot2","nr_spikes_u","nr_spikes_v","nr_spikes_w",
                     "nr_spikes_Ts","nr_spikes_h2o","nr_spikes_co2",
                     "ampl_res_u","ampl_res_v","ampl_res_w","ampl_res_Ts","ampl_res_h2o","ampl_res_co2",
                     "qf_most","qf_stationarity","qf_w","qf_all"]
    if do_methane:
        row["ch4_mean"] = ch4_mean
        row["ch4_sd"] = ch4_sd
        row["CH4F"] = cov_ch4w
        row["nr_spikes_ch4"] = nr_spikes_ch4
        row["ampl_res_ch4"] = ampl_res_ch4
    print(row)
    if do_log:
        logging.debug("Output row:\n%s", row)
    return(row)