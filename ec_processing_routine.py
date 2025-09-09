import numpy as np

def ec_processing_rt(out):
    """
    eddy-covariance post-processing of high-frequency data from Finse 

    @input numpy array with columns: u-v-w-Ts-h2o-co2 (high-frequency data from table finsefluxHFData)
    @return post-processed, quality-controlled, averaged fluxes and co (for inserting back to table ch_finsefluxpostproc)
    
    h2o, co2: used molar concentrations [mol/m^3] and converted them to density [kg/m^3] by multiplying with molar mass
    
    """
    #units
    out[:,3]=sos2Ts(out[:,3])
    #out[:,4]=ppt2rho(out[:,4],gas="H2O") #ppt
    #out[:,5]=ppt2rho(out[:,5]/1000,gas="CO2") #ppm originally
    out[:,4]=molarconcentration2density(out[:,4]/1000,gas="H2O") #mmol/m^3 -> kg/m^3
    out[:,5]=molarconcentration2density(out[:,5]/1000,gas="CO2")
    #despiking + spike counting
    nr_spikes_u=count_spikes(out[:,0],-15,15)
    out[:,0]=despiking(out[:,0],-15,15)
    nr_spikes_v=count_spikes(out[:,1],-15,15)
    out[:,1]=despiking(out[:,1],-15,15)
    nr_spikes_w=count_spikes(out[:,2],-5,5)
    out[:,2]=despiking(out[:,2],-5,5)
    nr_spikes_Ts=count_spikes(out[:,3],230,300)
    out[:,3]=despiking(out[:,3],230,300)
    nr_spikes_h2o=count_spikes(out[:,4],0,0.05)
    out[:,4]=despiking(out[:,4],0,0.05)
    nr_spikes_co2=count_spikes(out[:,5],0,0.005)
    out[:,5]=despiking(out[:,5],0,0.005)
    #amplitude resolution
    ampl_res_u=get_amplitude_resolution(out[:,0])
    ampl_res_v=get_amplitude_resolution(out[:,1])
    ampl_res_w=get_amplitude_resolution(out[:,2])
    ampl_res_Ts=get_amplitude_resolution(out[:,3])
    ampl_res_h2o=get_amplitude_resolution(out[:,4])
    ampl_res_co2=get_amplitude_resolution(out[:,5])
    #wind direction + mean
    wd_mean=calc_circular_mean(calc_windDirection(out[:,0],out[:,1]))
    ws_mean=np.nanmean(calc_windspeed(out[:,0],out[:,1]))
    #double rotation
    wind_rotated=rotate_double(out[:,0],out[:,1],out[:,2])
    out[:,0],out[:,1],out[:,2]=wind_rotated[0:3]
    dr_rot1=wind_rotated[3]
    dr_rot2=wind_rotated[4]
    #eddy stats use pandas quickly
    out=pd.DataFrame(out)
    #out=out.asfreq(freq='200L')
    u_mean,v_mean,w_mean,Ts_mean,h2o_mean,co2_mean=out.apply(np.mean)
    u_sd,v_sd,w_sd,Ts_sd,h2o_sd,co2_sd=out.apply(np.std)
    #back to numpy
    out=np.array(out)
    q_mean=np.nanmean(out[:,4])/1.225 #kg/kg
    T_mean=Ts2T(Ts_mean,q_mean)
    cov_uw=np.cov(out[:,0],out[:,2])[0,1]
    cov_uv=np.cov(out[:,0],out[:,1])[0,1]
    cov_vw=np.cov(out[:,1],out[:,2])[0,1]
    cov_wTs=np.cov(out[:,2],out[:,3])[0,1]
    cov_h2ow=np.cov(out[:,2],out[:,4])[0,1]
    cov_co2w=np.cov(out[:,2],out[:,5])[0,1]
    #snd correction
    cov_wT=SNDcorrection(Ts_mean,u_mean,v_mean,cov_uw,cov_vw,cov_wTs)
    #wpl correction
    #cov_h2ow=WPLcorrection(Ts_mean,q_mean,cov_wTs,h2o_mean,cov_h2ow)
    #cov_co2w=WPLcorrection(Ts_mean,q_mean,cov_wTs,h2o_mean,cov_h2ow,co2_mean,cov_co2w)
    #some more internal turbulence quantities
    sh=cov2sh(cov_wT)
    lh=cov2lh(cov_h2ow)
    et=lh2et(lh)
    br=calc_br(sh,lh)
    cf=cov2cf(cov_co2w)
    tke=calc_tke(u_sd,v_sd,w_sd)
    ustar=calc_ustar(cov_uw,cov_vw)
    z0=ustar2z0(ustar)
    dshear=calc_dshear(cov_uw,cov_vw)
    obukhov_length=calc_L(ustar,T_mean,cov_wT)
    zeta=calc_zeta(4.4,obukhov_length)
    #anisotropy todo
    xb=0
    yb=0
    flux_intermittency=0
    #quality flags
    qf_w=flag_w(w_mean)
    qf_most=flag_most(w_sd,ustar,zeta)
    #todo stationarity flag
    qf_stationarity=0
    qf_all=max(qf_w,qf_most,qf_stationarity)
    ### plots
    #qa
    #mrd
    #flux footprint
    #anisotropy
    ### return
    return(np.array([u_mean,v_mean,w_mean,Ts_mean,T_mean,h2o_mean,co2_mean,
                     u_sd,v_sd,w_sd,Ts_sd,h2o_sd,co2_sd,
                     cov_uw,cov_vw,cov_uv,wd_mean,ws_mean,
                     ustar,tke,dshear,z0,
                     sh,lh,et,br,cf,
                     obukhov_length,zeta,xb,yb,flux_intermittency,
                     np.int32(dr_rot1),np.int32(dr_rot2),nr_spikes_u,nr_spikes_v,nr_spikes_w,
                     nr_spikes_Ts,nr_spikes_h2o,nr_spikes_co2,
                     ampl_res_u,ampl_res_v,ampl_res_w,ampl_res_Ts,ampl_res_h2o,ampl_res_co2,
                     qf_most,qf_stationarity,qf_w,qf_all]))