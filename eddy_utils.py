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
    Vd=Runiversal*T_mean/(pres-e) #volume of dry air [m^3/mol]
    if gas == "H2O":
        return ppt/1000*M_H2O/Vd
    elif gas == "CO2":
        return ppt/1000*M_CO2/Vd
    elif gas == "CH4":
        return ppt/1000*M_CH4/Vd
